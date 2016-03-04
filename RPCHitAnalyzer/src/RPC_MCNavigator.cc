#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/RPC_MCNavigator.h"

namespace RPC
{
  std::vector<SimTrack>::const_iterator getTrack(const std::vector<SimTrack> &simTracks, unsigned int itrack)
  {
    for (std::vector<SimTrack>::const_iterator itt= simTracks.begin(); itt != simTracks.end(); itt++)
      if (itt->trackId() == itrack) return itt;
    return simTracks.end();
  }

  std::vector<SimVertex>::const_iterator getVertex(const std::vector<SimVertex> &simVertices, const SimTrack& track)
  {
    if (track.noVertex()) return simVertices.end();
    unsigned int ivertex=track.vertIndex();
    for (std::vector<SimVertex>::const_iterator itv= simVertices.begin(); itv != simVertices.end(); itv++)
      if (itv->vertexId() == ivertex) return itv;
    return simVertices.end();
  }

  const SimTrack& getParentTrack(const std::vector<SimTrack> &simTracks, const std::vector<SimVertex> &simVertices, const SimTrack& track)
  {
    if (! track.noGenpart()) return track;
    if (track.noVertex()) return track;
    std::vector<SimVertex>::const_iterator ivert=getVertex(simVertices,track);
    if (ivert == simVertices.end() ) return track;
    if (ivert->noParent()) return track;
    std::vector<SimTrack>::const_iterator itrack=getTrack(simTracks,ivert->parentIndex());
    if (itrack == simTracks.end()) return track;
    return *itrack;
  }

}
