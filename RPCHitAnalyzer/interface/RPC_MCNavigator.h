//helper function to navigate within the MC history
#ifndef RPC_MCNAVIGATOR_HH
#define RPC_MCNAVIGATOR_HH

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"

#include <vector>

namespace RPC
{

  //global  to define c speed (FIXME : might have that elsewhere in CMSSW)
  inline double cspeed() {return 29.9792458; /* cm/ns */}

  std::vector<SimTrack>::const_iterator getTrack(const std::vector<SimTrack> &simTracks, unsigned int trackid);

  inline std::vector<SimTrack>::const_iterator getTrack(const std::vector<SimTrack> &simTracks, const PSimHit& hit)
    {return getTrack(simTracks,hit.trackId());}

  std::vector<SimVertex>::const_iterator getVertex(const std::vector<SimVertex> &simVertices, const SimTrack& track);

  const SimTrack& getParentTrack(const std::vector<SimTrack> &simTracks, const std::vector<SimVertex> &simVertices, const SimTrack& track);
  
}


#endif
