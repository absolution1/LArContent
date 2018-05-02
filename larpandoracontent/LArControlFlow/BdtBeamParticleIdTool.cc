/**
 *  @file   larpandoracontent/LArControlFlow/BdtBeamParticleIdTool.cc
 *
 *  @brief  Implementation of the beam particle id tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/BdtBeamParticleIdTool.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

BdtBeamParticleIdTool::BdtBeamParticleIdTool() :
    m_useTrainingMode(false),
    m_trainingOutputFile(""),
    m_minPurity(0.8f),
    m_minCompleteness(0.8f),
    m_adaBoostDecisionTree(AdaBoostDecisionTree()),
    m_filePathEnvironmentVariable("FW_SEARCH_PATH"),
    m_maxNeutrinos(std::numeric_limits<int>::max()),
    m_minAdaBDTScore(0.f),
    m_sliceFeatureParameters(SliceFeatureParameters())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BdtBeamParticleIdTool::Initialize()
{
    // Get global LArTPC geometry information
    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);
    float parentMinX(pFirstLArTPC->GetCenterX() - 0.5f * pFirstLArTPC->GetWidthX());
    float parentMaxX(pFirstLArTPC->GetCenterX() + 0.5f * pFirstLArTPC->GetWidthX());
    float parentMinY(pFirstLArTPC->GetCenterY() - 0.5f * pFirstLArTPC->GetWidthY());
    float parentMaxY(pFirstLArTPC->GetCenterY() + 0.5f * pFirstLArTPC->GetWidthY());
    float parentMinZ(pFirstLArTPC->GetCenterZ() - 0.5f * pFirstLArTPC->GetWidthZ());
    float parentMaxZ(pFirstLArTPC->GetCenterZ() + 0.5f * pFirstLArTPC->GetWidthZ());

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);
        parentMinX = std::min(parentMinX, pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX());
        parentMaxX = std::max(parentMaxX, pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX());
        parentMinY = std::min(parentMinY, pLArTPC->GetCenterY() - 0.5f * pLArTPC->GetWidthY());
        parentMaxY = std::max(parentMaxY, pLArTPC->GetCenterY() + 0.5f * pLArTPC->GetWidthY());
        parentMinZ = std::min(parentMinZ, pLArTPC->GetCenterZ() - 0.5f * pLArTPC->GetWidthZ());
        parentMaxZ = std::max(parentMaxZ, pLArTPC->GetCenterZ() + 0.5f * pLArTPC->GetWidthZ());
    }

    try
    {
        m_sliceFeatureParameters.Initialize(parentMinX, parentMaxX, parentMinY, parentMaxY, parentMinZ, parentMaxZ);
    }
    catch (const StatusCodeException &statusCodeException)
    {
        std::cout << "BdtBeamParticleIdTool::Initialize - unable to initialize LArTPC geometry parameters" << std::endl;
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::SelectOutputPfos(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, PfoList &selectedPfos)
{
    if (nuSliceHypotheses.size() != crSliceHypotheses.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const unsigned int nSlices(nuSliceHypotheses.size());
    if (nSlices == 0) return;

    SliceFeaturesVector sliceFeaturesVector;
    this->GetSliceFeatures(this, nuSliceHypotheses, crSliceHypotheses, sliceFeaturesVector);

    if (m_useTrainingMode)
    {
        // ATTN in training mode, just return everything as a cosmic-ray
        this->SelectAllPfos(crSliceHypotheses, selectedPfos);

        pandora::IntVector bestSliceIndices;
        this->GetBestMCSliceIndices(pAlgorithm, nuSliceHypotheses, crSliceHypotheses, bestSliceIndices);

        for (unsigned int sliceIndex = 0; sliceIndex < nSlices; ++sliceIndex)
        {
            const SliceFeatures &features(sliceFeaturesVector.at(sliceIndex));
            if (!features.IsFeatureVectorAvailable()) continue;

            LArMvaHelper::MvaFeatureVector featureVector;

            try
            {
                features.FillFeatureVector(featureVector);
            }
            catch (const StatusCodeException &statusCodeException)
            {
                std::cout << "BdtBeamParticleIdTool::SelectOutputPfos - unable to fill feature vector" << std::endl;
                continue;
            }

            bool isGoodTrainingSlice(false);
            if (std::find(bestSliceIndices.begin(), bestSliceIndices.end(), sliceIndex) != bestSliceIndices.end())
                isGoodTrainingSlice = true;

            LArMvaHelper::ProduceTrainingExample(m_trainingOutputFile, isGoodTrainingSlice, featureVector);
        }

        return;
    }

    this->SelectPfosByAdaBDTScore(nuSliceHypotheses, crSliceHypotheses, sliceFeaturesVector, selectedPfos);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::GetSliceFeatures(const BdtBeamParticleIdTool *const pTool, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, SliceFeaturesVector &sliceFeaturesVector) const
{
    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
        sliceFeaturesVector.push_back(SliceFeatures(nuSliceHypotheses.at(sliceIndex), crSliceHypotheses.at(sliceIndex), pTool, m_sliceFeatureParameters));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::SelectAllPfos(const SliceHypotheses &hypotheses, PfoList &selectedPfos) const
{
    for (const PfoList &pfos : hypotheses)
        this->SelectPfos(pfos, selectedPfos);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::SelectPfos(const PfoList &pfos, PfoList &selectedPfos) const
{
    selectedPfos.insert(selectedPfos.end(), pfos.begin(), pfos.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::GetBestMCSliceIndices(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, pandora::IntVector &bestSliceIndices) const
{
    // Get all hits in all slices to find true number of mc hits
    const CaloHitList *pAllReconstructedCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, m_caloHitListName, pAllReconstructedCaloHitList));

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, m_mcParticleListName, pMCParticleList));

    // Obtain map: [mc particle -> primary mc particle]
    LArMCParticleHelper::MCRelationMap mcToPrimaryMCMap;
    LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToPrimaryMCMap);

    // Remove non-reconstructable hits, e.g. those downstream of a neutron
    CaloHitList reconstructableCaloHitList;
    LArMCParticleHelper::PrimaryParameters parameters;
    LArMCParticleHelper::SelectCaloHits(pAllReconstructedCaloHitList, mcToPrimaryMCMap, reconstructableCaloHitList, parameters.m_selectInputHits, parameters.m_maxPhotonPropagation);

    MCParticleToIntMap mcParticleToReconstructableHitsMap;
    this->PopulateMCParticleToHitsMap(mcParticleToReconstructableHitsMap, reconstructableCaloHitList);

    const CaloHitSet reconstructableCaloHitSet(reconstructableCaloHitList.begin(), reconstructableCaloHitList.end());

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        // All hits in a slice - No double counting
        CaloHitList reconstructedCaloHitList;
        this->Collect2DHits(crSliceHypotheses.at(sliceIndex), reconstructedCaloHitList, reconstructableCaloHitSet);

        if (nuSliceHypotheses.at(sliceIndex).size() == 1)
        {
            const PfoList &nuFinalStates(nuSliceHypotheses.at(sliceIndex).front()->GetDaughterPfoList());
            this->Collect2DHits(nuFinalStates, reconstructedCaloHitList, reconstructableCaloHitSet);
        }

        const unsigned int nRecoHits(reconstructedCaloHitList.size());

        // MCParticle to hits in slice map
        MCParticleToIntMap mcParticleToHitsInSliceMap;
        this->PopulateMCParticleToHitsMap(mcParticleToHitsInSliceMap, reconstructedCaloHitList);

        if (mcParticleToHitsInSliceMap.empty())
            continue;

        // Get best mc particle for slice
        const MCParticle *pBestMCParticle(nullptr);
        unsigned int nSharedHits(0);

        for (const auto &iter : mcParticleToHitsInSliceMap)
        {
            if (iter.second > static_cast<int>(nSharedHits))
            {
                pBestMCParticle = iter.first;
                nSharedHits = iter.second;
            }
        }

        // Only consider if target beam particles
        if (2001 != LArMCParticleHelper::GetNuanceCode(pBestMCParticle))
            continue;

        const unsigned int nMCHits(mcParticleToReconstructableHitsMap.at(pBestMCParticle));
        const float purity(nRecoHits > 0 ? static_cast<float>(nSharedHits) / static_cast<float>(nRecoHits) : 0.f);
        const float completeness(nMCHits > 0 ? static_cast<float>(nSharedHits) / static_cast<float>(nMCHits) : 0.f);

        if (this->PassesQualityCuts(purity, completeness))
            bestSliceIndices.push_back(sliceIndex);
    }
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::PopulateMCParticleToHitsMap(MCParticleToIntMap &mcParticleToIntMap, const pandora::CaloHitList &caloHitList) const
{
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        try
        {
            const MCParticle *pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit)));
            MCParticleToIntMap::iterator iter(mcParticleToIntMap.find(pParentMCParticle));

            if (iter != mcParticleToIntMap.end())
            {
                (*iter).second += 1;
            }
            else
            {
                mcParticleToIntMap.insert(MCParticleToIntMap::value_type(pParentMCParticle, 1));
            }
        }
        catch (const StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_NOT_INITIALIZED != statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::Collect2DHits(const PfoList &pfos, CaloHitList &reconstructedCaloHitList, const CaloHitSet &reconstructableCaloHitSet) const
{
    CaloHitList collectedHits;
    LArPfoHelper::GetCaloHits(pfos, TPC_VIEW_U, collectedHits);
    LArPfoHelper::GetCaloHits(pfos, TPC_VIEW_V, collectedHits);
    LArPfoHelper::GetCaloHits(pfos, TPC_VIEW_W, collectedHits);

    for (const CaloHit *const pCaloHit : collectedHits)
    {
        // ATTN hits collected from Pfos are copies of hits passed from master instance, we need to access their parent to use MC info
        const CaloHit *pParentCaloHit(static_cast<const CaloHit *>(pCaloHit->GetParentAddress()));

        if (!reconstructableCaloHitSet.count(pParentCaloHit))
            continue;

        // Ensure no hits have been double counted
        if (std::find(reconstructedCaloHitList.begin(), reconstructedCaloHitList.end(), pParentCaloHit) == reconstructedCaloHitList.end())
            reconstructedCaloHitList.push_back(pParentCaloHit); 
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool BdtBeamParticleIdTool::PassesQualityCuts(const float purity, const float completeness) const
{
    if ((purity < m_minPurity) || (completeness < m_minCompleteness)) return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::SelectPfosByAdaBDTScore(const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, const SliceFeaturesVector &sliceFeaturesVector, PfoList &selectedPfos) const
{
    // Calculate the probability of each slice that passes the minimum probability cut
    std::vector<UintFloatPair> sliceIndexAdaBDTScorePairs;
    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        const float nuAdaBDTScore(sliceFeaturesVector.at(sliceIndex).GetAdaBoostDecisionTreeScore(m_adaBoostDecisionTree));

        if (nuAdaBDTScore < m_minAdaBDTScore)
        {
            this->SelectPfos(crSliceHypotheses.at(sliceIndex), selectedPfos);
            continue;
        }

        sliceIndexAdaBDTScorePairs.push_back(UintFloatPair(sliceIndex, nuAdaBDTScore));
    }

    // Sort the slices by probability
    std::sort(sliceIndexAdaBDTScorePairs.begin(), sliceIndexAdaBDTScorePairs.end(), [] (const UintFloatPair &a, const UintFloatPair &b)
    {
        return (a.second > b.second);
    });

    // Select the first m_maxNeutrinos as neutrinos, and the rest as cosmic
    unsigned int nNuSlices(0);
    for (const UintFloatPair &slice : sliceIndexAdaBDTScorePairs)
    {
        if (nNuSlices < m_maxNeutrinos)
        {
            this->SelectPfos(nuSliceHypotheses.at(slice.first), selectedPfos);
            nNuSlices++;
            continue;
        }

        this->SelectPfos(crSliceHypotheses.at(slice.first), selectedPfos);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

BdtBeamParticleIdTool::Plane::Plane(const CartesianVector &normal, const CartesianVector &point) :
    m_unitNormal(0.f, 0.f, 0.f),
    m_point(point),
    m_d(-1.f * (normal.GetDotProduct(point)))
{
    try
    {
        m_unitNormal = normal.GetUnitVector();
    }
    catch (const StatusCodeException &statusCodeException)
    {
        std::cout << "BdtBeamParticleIdTool::Plane::Plane - normal vector to plane has a magnitude of zero" << std::endl;
        throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector BdtBeamParticleIdTool::Plane::GetLineIntersection(const CartesianVector &point, const CartesianVector &direction) const
{
    if (std::fabs(direction.GetDotProduct(m_unitNormal)) < std::numeric_limits<float>::epsilon())
        return CartesianVector(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());

    const float denominator(direction.GetDotProduct(m_unitNormal));

    if (std::fabs(denominator) < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

    const float t(-1.f * (point.GetDotProduct(m_unitNormal) + m_d) / denominator);
    return (point + (direction * t));
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

BdtBeamParticleIdTool::SliceFeatureParameters::SliceFeatureParameters() :
    m_beamLArTPCIntersection(0.f, 0.f, 0.f),
    m_beamDirection(0.f, 0.f, 0.f),
    m_selectedFraction(10.f),
    m_nSelectedHits(100)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::SliceFeatureParameters::Initialize(const float larTPCMinX, const float larTPCMaxX, const float larTPCMinY, const float larTPCMaxY, const float larTPCMinZ, const float larTPCMaxZ)
{
    m_larTPCMinX = larTPCMinX;
    m_larTPCMaxX = larTPCMaxX;
    m_larTPCMinY = larTPCMinY;
    m_larTPCMaxY = larTPCMaxY;
    m_larTPCMinZ = larTPCMinZ;
    m_larTPCMaxZ = larTPCMaxZ;

    const CartesianVector normalTop(0.f, 0.f, 1.f), pointTop(0.f, 0.f, m_larTPCMaxZ);
    const CartesianVector normalBottom(0.f, 0.f, -1.f), pointBottom(0.f, 0.f, m_larTPCMinZ);
    const CartesianVector normalRight(1.f, 0.f, 0.f), pointRight(m_larTPCMaxX, 0.f, 0.f);
    const CartesianVector normalLeft(-1.f, 0.f, 0.f), pointLeft(m_larTPCMinX, 0.f, 0.f);
    const CartesianVector normalBack(0.f, 1.f, 0.f), pointBack(0.f, m_larTPCMaxY, 0.f);
    const CartesianVector normalFront(0.f, -1.f, 0.f), pointFront(0.f, m_larTPCMinY, 0.f);

    m_larTPCPlanes.emplace_back(normalTop, pointTop);
    m_larTPCPlanes.emplace_back(normalBottom, pointBottom);
    m_larTPCPlanes.emplace_back(normalRight, pointRight);
    m_larTPCPlanes.emplace_back(normalLeft, pointLeft);
    m_larTPCPlanes.emplace_back(normalBack, pointBack);
    m_larTPCPlanes.emplace_back(normalFront, pointFront);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

BdtBeamParticleIdTool::SliceFeatures::SliceFeatures(const PfoList &pfosNu, const PfoList &pfosCr, const BdtBeamParticleIdTool *const pTool, const SliceFeatureParameters &sliceFeatureParameters) :
    m_isAvailable(false),
    m_sliceFeatureParameters(sliceFeatureParameters),
    m_pTool(pTool)
{
    try
    {
        double closestDistanceNu(std::numeric_limits<double>::max()), closestDistanceCr(std::numeric_limits<double>::max()), supplementaryAngleToBeamNu(std::numeric_limits<double>::max()),
               supplementaryAngleToBeamCr(std::numeric_limits<double>::max()), separationNu(std::numeric_limits<double>::max()), separationCr(std::numeric_limits<double>::max());;
        CaloHitList caloHitList3DNu, caloHitList3DCr, selectedCaloHitListNu, selectedCaloHitListCr;
        PfoList allConnectedPfoListNu, allConnectedPfoListCr;
        LArPcaHelper::EigenValues eigenValuesNu(0.f, 0.f, 0.f);
        LArPcaHelper::EigenValues eigenValuesCr(0.f, 0.f, 0.f);

        LArPfoHelper::GetAllConnectedPfos(pfosNu, allConnectedPfoListNu);
        LArPfoHelper::GetAllConnectedPfos(pfosCr, allConnectedPfoListCr);
        LArPfoHelper::GetCaloHits(allConnectedPfoListNu, TPC_3D, caloHitList3DNu);
        LArPfoHelper::GetCaloHits(allConnectedPfoListCr, TPC_3D, caloHitList3DCr);

        this->GetLeadingCaloHits(caloHitList3DNu, selectedCaloHitListNu, closestDistanceNu);
        this->GetLeadingCaloHits(caloHitList3DCr, selectedCaloHitListCr, closestDistanceCr);

        if (!selectedCaloHitListNu.empty() && !selectedCaloHitListCr.empty())
        {
            float maxYNu(-std::numeric_limits<float>::max()), maxYCr(-std::numeric_limits<float>::max());
            CartesianVector centroidNu(0.f, 0.f, 0.f), interceptOneNu(0.f, 0.f, 0.f), interceptTwoNu(0.f, 0.f, 0.f), centroidCr(0.f, 0.f, 0.f), interceptOneCr(0.f, 0.f, 0.f), interceptTwoCr(0.f, 0.f, 0.f);

            // Beam
            LArPcaHelper::EigenVectors eigenVecsNu;
            LArPcaHelper::RunPca(selectedCaloHitListNu, centroidNu, eigenValuesNu, eigenVecsNu);
            const CartesianVector &majorAxisNu(eigenVecsNu.front());
            supplementaryAngleToBeamNu = majorAxisNu.GetOpeningAngle(m_sliceFeatureParameters.GetBeamDirection());

            this->GetLArTPCIntercepts(centroidNu, majorAxisNu, interceptOneNu, interceptTwoNu);
            const double separationOneNu((interceptOneNu - m_sliceFeatureParameters.GetBeamLArTPCIntersection()).GetMagnitude());
            const double separationTwoNu((interceptTwoNu - m_sliceFeatureParameters.GetBeamLArTPCIntersection()).GetMagnitude());
            separationNu = std::min(separationOneNu, separationTwoNu);

            // Cosmic
            LArPcaHelper::EigenVectors eigenVecsCr;
            LArPcaHelper::RunPca(selectedCaloHitListCr, centroidCr, eigenValuesCr, eigenVecsCr);
            const CartesianVector &majorAxisCr(eigenVecsCr.front());
            supplementaryAngleToBeamCr = majorAxisCr.GetOpeningAngle(m_sliceFeatureParameters.GetBeamDirection());

            this->GetLArTPCIntercepts(centroidCr, majorAxisCr, interceptOneCr, interceptTwoCr);
            const double separationOneCr((interceptOneCr - m_sliceFeatureParameters.GetBeamLArTPCIntersection()).GetMagnitude());
            const double separationTwoCr((interceptTwoCr - m_sliceFeatureParameters.GetBeamLArTPCIntersection()).GetMagnitude());
            separationCr = std::min(separationOneCr, separationTwoCr);

            for (const CaloHit *pCaloHit : caloHitList3DNu)
            {
                if (maxYNu < pCaloHit->GetPositionVector().GetY())
                    maxYNu = pCaloHit->GetPositionVector().GetY();
            }

            for (const CaloHit *pCaloHit : caloHitList3DCr)
            {
                if (maxYCr < pCaloHit->GetPositionVector().GetY())
                    maxYCr = pCaloHit->GetPositionVector().GetY();
            }

            m_featureVector.push_back(closestDistanceNu);
            m_featureVector.push_back(supplementaryAngleToBeamNu);
            m_featureVector.push_back(separationNu);
            m_featureVector.push_back(eigenValuesNu.GetX());
            m_featureVector.push_back(eigenValuesNu.GetY());
            m_featureVector.push_back(eigenValuesNu.GetZ());
            m_featureVector.push_back(maxYNu);
            m_featureVector.push_back(allConnectedPfoListNu.size());

            m_featureVector.push_back(closestDistanceCr);
            m_featureVector.push_back(supplementaryAngleToBeamCr);
            m_featureVector.push_back(separationCr);
            m_featureVector.push_back(eigenValuesCr.GetX());
            m_featureVector.push_back(eigenValuesCr.GetY());
            m_featureVector.push_back(eigenValuesCr.GetZ());
            m_featureVector.push_back(maxYCr);
            m_featureVector.push_back(allConnectedPfoListCr.size());

            m_isAvailable = true;
        }
    }
    catch (const StatusCodeException &)
    {
        return;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::SliceFeatures::GetLeadingCaloHits(const CaloHitList &inputCaloHitList, CaloHitList &outputCaloHitList,
    double &closestHitToFaceDistance) const
{
    if (inputCaloHitList.empty())
    {
        std::cout << "BdtBeamParticleIdTool::SliceFeatures::GetLeadingCaloHits - empty calo hit list" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
    }

    typedef std::pair<const CaloHit*, float> HitDistancePair;
    typedef std::vector<HitDistancePair> HitDistanceVector;
    HitDistanceVector hitDistanceVector;

    for (const CaloHit *const pCaloHit : inputCaloHitList)
        hitDistanceVector.emplace_back(pCaloHit, (pCaloHit->GetPositionVector() - m_sliceFeatureParameters.GetBeamLArTPCIntersection()).GetMagnitudeSquared());

    std::sort(hitDistanceVector.begin(), hitDistanceVector.end(), [](const HitDistancePair &lhs, const HitDistancePair &rhs) -> bool {return (lhs.second < rhs.second);});

    if (hitDistanceVector.front().second < 0.f)
    {
        std::cout << "BdtBeamParticleIdTool::SliceFeatures::GetLeadingCaloHits - unphysical magnitude of a vector" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
    }

    closestHitToFaceDistance = std::sqrt(hitDistanceVector.front().second);

    const unsigned int nInputHits(inputCaloHitList.size());
    const unsigned int nSelectedCaloHits(nInputHits < m_sliceFeatureParameters.GetNSelectedHits() ? nInputHits :
        static_cast<unsigned int>(std::ceil(static_cast<float>(nInputHits) * m_sliceFeatureParameters.GetSelectedFraction() / 100.f)));

    for (const HitDistancePair &hitDistancePair : hitDistanceVector)
    {
        outputCaloHitList.push_back(hitDistancePair.first);

        if (outputCaloHitList.size() >= nSelectedCaloHits)
            break;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::SliceFeatures::GetLArTPCIntercepts(const CartesianVector &a0, const CartesianVector &lineDirection,
    CartesianVector &interceptOne, CartesianVector &interceptTwo) const
{
    CartesianPointVector intercepts;
    CartesianVector lineUnitVector(0.f, 0.f, 0.f);

    try
    {
        lineUnitVector = lineDirection.GetUnitVector();
    }
    catch (const StatusCodeException &statusCodeException)
    {
        std::cout << "BdtBeamParticleIdTool::SliceFeatures::GetLArTPCIntercepts - normal vector to plane has a magnitude of zero" << std::endl;
        throw statusCodeException;
    }

    for (const Plane &plane : m_sliceFeatureParameters.GetPlanes())
    {
        const CartesianVector intercept(plane.GetLineIntersection(a0, lineUnitVector));

        if (this->IsContained(intercept))
            intercepts.push_back(intercept);
    }

    if (intercepts.size() == 2)
    {
        interceptOne = intercepts.at(0);
        interceptTwo = intercepts.at(1);
    }
    else
    {
        std::cout << "BdtBeamParticleIdTool::SliceFeatures::GetLArTPCIntercepts - inconsistent number of intercepts between a line and the LArTPC" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool BdtBeamParticleIdTool::SliceFeatures::IsContained(const CartesianVector &spacePoint) const
{
    if ((m_sliceFeatureParameters.GetLArTPCMinX() - spacePoint.GetX() > std::numeric_limits<float>::epsilon()) || (spacePoint.GetX() - m_sliceFeatureParameters.GetLArTPCMaxX() > std::numeric_limits<float>::epsilon()) ||
        (m_sliceFeatureParameters.GetLArTPCMinY() - spacePoint.GetY() > std::numeric_limits<float>::epsilon()) || (spacePoint.GetY() - m_sliceFeatureParameters.GetLArTPCMaxY() > std::numeric_limits<float>::epsilon()) ||
        (m_sliceFeatureParameters.GetLArTPCMinZ() - spacePoint.GetZ() > std::numeric_limits<float>::epsilon()) || (spacePoint.GetZ() - m_sliceFeatureParameters.GetLArTPCMaxZ() > std::numeric_limits<float>::epsilon()))
    {
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::SliceFeatures::FillFeatureVector(LArMvaHelper::MvaFeatureVector &featureVector) const
{
    if (!m_isAvailable)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    if (!featureVector.empty())
    {
        std::cout << "BdtBeamParticleIdTool::SliceFeatures::FillFeatureVector - feature vector already populated" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
    }

    featureVector.insert(featureVector.end(), m_featureVector.begin(), m_featureVector.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

float BdtBeamParticleIdTool::SliceFeatures::GetAdaBoostDecisionTreeScore(const AdaBoostDecisionTree &adaBoostDecisionTree) const
{
    // ATTN if one or more of the features can not be calculated, then default to calling the slice a cosmic ray.  -1.f is the minimum score 
    // possible for a weighted bdt.
    if (!this->IsFeatureVectorAvailable()) return -1.f;

    LArMvaHelper::MvaFeatureVector featureVector;

    try
    {
        this->FillFeatureVector(featureVector);
    }
    catch (const StatusCodeException &statusCodeException)
    {
        std::cout << "BdtBeamParticleIdTool::SliceFeatures::GetAdaBoostDecisionTreeScore - unable to fill feature vector" << std::endl;
        return -1.f;
    }

    return LArMvaHelper::CalculateClassificationScore(adaBoostDecisionTree, featureVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BdtBeamParticleIdTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    // BDT Settings
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseTrainingMode", m_useTrainingMode));

    if (m_useTrainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "TrainingOutputFileName", m_trainingOutputFile));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "CaloHitListName", m_caloHitListName));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "MCParticleListName", m_mcParticleListName));
    }
    else
    {
        std::string bdtName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "BdtName", bdtName));

        std::string bdtFileName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "BdtFileName", bdtFileName));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "MinAdaBDTScore", m_minAdaBDTScore));

        const std::string fullBdtFileName(LArFileHelper::FindFileInPath(bdtFileName, m_filePathEnvironmentVariable));
        const StatusCode statusCode(m_adaBoostDecisionTree.Initialize(fullBdtFileName, bdtName));

        if (STATUS_CODE_SUCCESS != statusCode)
        {
            std::cout << "BdtBeamParticleIdTool::ReadSettings - unable to load bdt" << std::endl;
            return statusCode;
        }
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinimumPurity", m_minPurity));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinimumCompleteness", m_minCompleteness));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FilePathEnvironmentVariable", m_filePathEnvironmentVariable));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaximumNeutrinos", m_maxNeutrinos));

    // Geometry Information for training
    FloatVector beamLArTPCIntersection;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "BeamTPCIntersection", beamLArTPCIntersection));

    if (3 == beamLArTPCIntersection.size())
    {
        pandora::CartesianVector beamLArTPCIntersectionCartesianVector(beamLArTPCIntersection.at(0), beamLArTPCIntersection.at(1), beamLArTPCIntersection.at(2));
        m_sliceFeatureParameters.SetBeamLArTPCIntersection(beamLArTPCIntersectionCartesianVector);
    }
    else if (!beamLArTPCIntersection.empty())
    {
        std::cout << "BdtBeamParticleIdTool::ReadSettings - invalid BeamTPCIntersection specified " << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }
    else
    {
        // Default for protoDUNE.
        pandora::CartesianVector beamLArTPCIntersectionCartesianVector(-33.051f, 461.06f, 0.f);
        m_sliceFeatureParameters.SetBeamLArTPCIntersection(beamLArTPCIntersectionCartesianVector);
    }

    FloatVector beamDirection;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "BeamDirection", beamDirection));

    if (3 == beamDirection.size())
    {
        CartesianVector beamDirectionCartesianVector(beamDirection.at(0), beamDirection.at(1), beamDirection.at(2));
        m_sliceFeatureParameters.SetBeamDirection(beamDirectionCartesianVector);
    }
    else if (!beamDirection.empty())
    {
        std::cout << "BdtBeamParticleIdTool::ReadSettings - invalid BeamDirection specified " << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }
    else
    {
        // Default for protoDUNE.
        const float thetaXZ0(-11.844f * M_PI / 180.f);
        CartesianVector beamDirectionCartesianVector(std::sin(thetaXZ0), 0.f, std::cos(thetaXZ0));
        m_sliceFeatureParameters.SetBeamDirection(beamDirectionCartesianVector);
    }

    float selectedFraction(0.f);

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectedFraction", selectedFraction));

    if (selectedFraction > std::numeric_limits<float>::epsilon())
        m_sliceFeatureParameters.SetSelectedFraction(selectedFraction);

    unsigned int nSelectedHits(0);

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NSelectedHits", nSelectedHits));

    if (nSelectedHits > 0)
        m_sliceFeatureParameters.SetNSelectedHits(nSelectedHits);

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content