/**
 *  @file   larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewTransverseTracksAlgorithm.cc
 *
 *  @brief  Implementation of the two view transverse tracks algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArObjects/LArDiscreteCumulativeDistribution.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArDiscreteCumulativeDistributionHelper.h"


#include "larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewTransverseTracksAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TwoViewTransverseTracksAlgorithm::TwoViewTransverseTracksAlgorithm() :
    m_nMaxMatrixToolRepeats(1000)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewTransverseTracksAlgorithm::CalculateOverlapResult(const Cluster *const pCluster1, const Cluster *const pCluster2)
{
    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    LArClusterHelper::GetClusterSpanX(pCluster1, xMin1, xMax1);
    LArClusterHelper::GetClusterSpanX(pCluster2, xMin2, xMax2);

    const float xSpan1(xMax1 - xMin1), xSpan2(xMax2 - xMin2);

    if ((xSpan1 < std::numeric_limits<float>::epsilon()) || (xSpan2 < std::numeric_limits<float>::epsilon()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const float xOverlapMin(std::max(xMin1, xMin2));
    const float xOverlapMax(std::min(xMax1, xMax2));
    const float xOverlap(xOverlapMax - xOverlapMin);
    TwoViewXOverlap twoViewXOverlap(xMin1, xMax1, xMin2, xMax2, xOverlap);

    float zMin1(0.f), zMax1(0.f);
    float zMin2(0.f), zMax2(0.f);
    LArClusterHelper::GetClusterSpanZ(pCluster1, xMin1, xMax1, zMin1, zMax1);
    LArClusterHelper::GetClusterSpanZ(pCluster2, xMin2, xMax2, zMin2, zMax2);
    CartesianVector boundingBoxMin1(xOverlapMin, 0.f, zMin1), boundingBoxMax1(xOverlapMax, 0.f, zMax1);
    CartesianVector boundingBoxMin2(xOverlapMin, 0.f, zMin2), boundingBoxMax2(xOverlapMax, 0.f, zMax2);
    pandora::CaloHitList overlapHits1, overlapHits2;
    LArClusterHelper::GetCaloHitListInBoundingBox(pCluster1, boundingBoxMin1, boundingBoxMax1, overlapHits1);
    LArClusterHelper::GetCaloHitListInBoundingBox(pCluster2, boundingBoxMin2, boundingBoxMax2, overlapHits2);

    DiscreteCumulativeDistribution disCumulDist1, disCumulDist2;
    LArDiscreteCumulativeDistributionHelper::CreateDistributionFromCaloHits(overlapHits1, disCumulDist1);
    LArDiscreteCumulativeDistributionHelper::CreateDistributionFromCaloHits(overlapHits2, disCumulDist2);
    float matchingScore(LArDiscreteCumulativeDistributionHelper::CalculatePValueWithKSTestStatistic(disCumulDist1, disCumulDist2));

    TwoViewTransverseOverlapResult twoViewTransverseOverlapResult(matchingScore, twoViewXOverlap);

    if (xOverlap > std::numeric_limits<float>::epsilon())
        m_overlapMatrix.SetOverlapResult(pCluster1, pCluster2, twoViewTransverseOverlapResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewTransverseTracksAlgorithm::ExamineOverlapContainer()
{
    unsigned int repeatCounter(0);

    for (MatrixToolVector::const_iterator iter = m_algorithmToolVector.begin(), iterEnd = m_algorithmToolVector.end(); iter != iterEnd; )
    {
        if ((*iter)->Run(this, m_overlapMatrix))
        {
            iter = m_algorithmToolVector.begin();

            if (++repeatCounter > m_nMaxMatrixToolRepeats)
                break;
        }
        else
        {
            ++iter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewTransverseTracksAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "TrackTools", algorithmToolVector));

    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        TransverseMatrixTool *const pTransverseMatrixTool(dynamic_cast<TransverseMatrixTool*>(*iter));

        if (!pTransverseMatrixTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pTransverseMatrixTool);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NMaxMatrixToolRepeats", m_nMaxMatrixToolRepeats));

    return TwoViewTrackMatchingAlgorithm<TwoViewTransverseOverlapResult>::ReadSettings(xmlHandle);
}

} // namespace lar_content
