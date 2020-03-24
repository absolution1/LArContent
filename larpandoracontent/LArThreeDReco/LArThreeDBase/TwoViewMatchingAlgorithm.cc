/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/TwoViewMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the matching base algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArObjects/LArTrackTwoViewOverlapResult.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/TwoViewMatchingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

template <typename T>
TwoViewMatchingAlgorithm<T>::TwoViewMatchingAlgorithm() :
    m_pInputClusterList1(nullptr),
    m_pInputClusterList2(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
TwoViewMatchingAlgorithm<T>::~TwoViewMatchingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void TwoViewMatchingAlgorithm<T>::UpdateForNewCluster(const Cluster *const pNewCluster)
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(pNewCluster));
    HitTypeToIndexMap::const_iterator iter = m_hitTypeToIndexMap.find(hitType);

    if ((m_hitTypeToIndexMap.end() == iter) || ((1 != iter->second) && (2 != iter->second)))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    ClusterList &clusterList((1 == iter->second) ? m_clusterList1 : m_clusterList2);

    if (clusterList.end() != std::find(clusterList.begin(), clusterList.end(), pNewCluster))
        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

    clusterList.push_back(pNewCluster);

    const ClusterList &clusterList2((1 == iter->second) ? m_clusterList2 : m_clusterList1);

    ClusterVector clusterVector2(clusterList2.begin(), clusterList2.end());
    std::sort(clusterVector2.begin(), clusterVector2.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster2 : clusterVector2)
    {
        if (1 == iter->second)
        {
            this->CalculateOverlapResult(pNewCluster, pCluster2);
        }
        else
        {
            this->CalculateOverlapResult(pCluster2, pNewCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void TwoViewMatchingAlgorithm<T>::UpdateUponDeletion(const Cluster *const pDeletedCluster)
{
    ClusterList::iterator iter1 = std::find(m_clusterList1.begin(), m_clusterList1.end(), pDeletedCluster);
    ClusterList::iterator iter2 = std::find(m_clusterList2.begin(), m_clusterList2.end(), pDeletedCluster);

    if (m_clusterList1.end() != iter1)
        m_clusterList1.erase(iter1);

    if (m_clusterList2.end() != iter2)
        m_clusterList2.erase(iter2);

    m_overlapMatrix.RemoveCluster(pDeletedCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
const std::string &TwoViewMatchingAlgorithm<T>::GetClusterListName(const HitType hitType) const
{
    HitTypeToIndexMap::const_iterator iter = m_hitTypeToIndexMap.find(hitType);
    if (m_hitTypeToIndexMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if ((1 != iter->second) && (2 != iter->second))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return ((1 == iter->second) ? m_inputClusterListName1 : m_inputClusterListName2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
const pandora::ClusterList &TwoViewMatchingAlgorithm<T>::GetInputClusterList(const HitType hitType) const
{
    HitTypeToIndexMap::const_iterator iter = m_hitTypeToIndexMap.find(hitType);
    if (m_hitTypeToIndexMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    if ((1 == iter->second) && m_pInputClusterList1)
        return (*m_pInputClusterList1);

    if ((2 == iter->second) && m_pInputClusterList2)
        return (*m_pInputClusterList2);

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
const pandora::ClusterList &TwoViewMatchingAlgorithm<T>::GetSelectedClusterList(const HitType hitType) const
{
    HitTypeToIndexMap::const_iterator iter = m_hitTypeToIndexMap.find(hitType);
    if (m_hitTypeToIndexMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if ((1 != iter->second) && (2 != iter->second))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return ((1 == iter->second) ? m_clusterList1 : m_clusterList2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void TwoViewMatchingAlgorithm<T>::SelectAllInputClusters()
{
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        m_inputClusterListName1, m_pInputClusterList1));
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        m_inputClusterListName2, m_pInputClusterList2));

    if (!m_pInputClusterList1 || !m_pInputClusterList2)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "TwoViewMatchingAlgorithm: one or more input cluster lists unavailable." << std::endl;

        throw StatusCodeException(STATUS_CODE_SUCCESS);
    }

    if (!m_pInputClusterList1->empty())
        m_hitTypeToIndexMap.insert(HitTypeToIndexMap::value_type(LArClusterHelper::GetClusterHitType(m_pInputClusterList1->front()), 1));

    if (!m_pInputClusterList2->empty())
        m_hitTypeToIndexMap.insert(HitTypeToIndexMap::value_type(LArClusterHelper::GetClusterHitType(m_pInputClusterList2->front()), 2));

    this->SelectInputClusters(m_pInputClusterList1, m_clusterList1);
    this->SelectInputClusters(m_pInputClusterList2, m_clusterList2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void TwoViewMatchingAlgorithm<T>::TidyUp()
{
    m_overlapMatrix.Clear();
    m_hitTypeToIndexMap.clear();

    m_pInputClusterList1 = nullptr;
    m_pInputClusterList2 = nullptr;

    m_clusterList1.clear();
    m_clusterList2.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void TwoViewMatchingAlgorithm<T>::PerformMainLoop()
{
    // TODO Support looping over all pairs of lists from an input vector of list names
    ClusterVector clusterVector1(m_clusterList1.begin(), m_clusterList1.end());
    ClusterVector clusterVector2(m_clusterList2.begin(), m_clusterList2.end());
    std::sort(clusterVector1.begin(), clusterVector1.end(), LArClusterHelper::SortByNHits);
    std::sort(clusterVector2.begin(), clusterVector2.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster1 : clusterVector1)
    {
        for (const Cluster *const pCluster2 : clusterVector2)
            this->CalculateOverlapResult(pCluster1, pCluster2);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
StatusCode TwoViewMatchingAlgorithm<T>::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListName1", m_inputClusterListName1));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListName2", m_inputClusterListName2));

    return MatchingBaseAlgorithm::ReadSettings(xmlHandle);
}

template class TwoViewMatchingAlgorithm<float>;
template class TwoViewMatchingAlgorithm<TwoViewTransverseOverlapResult>;

} // namespace lar_content
