#ifndef IBDALG_HH
#define IBDALG_HH

#include "SniperKernel/AlgBase.h"
#include <vector>
#include <algorithm>
#include <unordered_map> 
#include <math.h>
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TVector3.h"
#include "TTimeStamp.h"

#include "EvtNavigator/NavBuffer.h"
#include "EvtNavigator/EvtNavigator.h"
#include "EvtNavigator/EvtNavHelper.h"
#include "Event/OecHeader.h"
#include "Event/OecEvt.h"
#include "OECTagSvc/OECTagSvc.h"
#include "OECTagID/OECTagID.h"

class IBDSelectionAlg : public AlgBase
{
private:
    JM::NavBuffer* m_buf;
    OECTagSvc* m_tagsvc;

    JM::NavBuffer::Iterator DelayIt;
    // JM::EvtNavigator* PostDelayNav;

    bool foundDelay;
    int m_DelayIBD;
    std::string m_DelayType;

    bool isMuonVetoed(JM::NavBuffer::Iterator);
    bool isIsolated(JM::NavBuffer::Iterator);
    bool findCorrelation(JM::NavBuffer::Iterator, std::vector<JM::OecEvt*>&);

public:
    IBDSelectionAlg(const std::string&);

    bool initialize();
    bool execute();
    bool finalize();

    bool bookTree();

private:

    TGraph* gInterfaceLevel;

    int m_iEvt;
    int m_iRun;

    uint32_t i_pBiPo214;
    uint32_t i_dBiPo214;

    unsigned int TotalLPMT = 17612;
    unsigned int TotalSPMT = 25600;
    
    std::vector<TVector3> ALL_LPMT_pos;
    std::vector<TVector3> ALL_SPMT_pos;

    uint64_t m_TimeStamp;

    std::vector<int> m_PmtId;
    std::vector<double> m_HitTime;
    std::vector<std::string> m_TriggerType;

    std::vector<double> m_PromptHitTimeTOF;
    std::vector<double> m_DelayHitTimeTOF;
    std::vector<std::string> m_EventTag;

    uint64_t m_TimeDifference;

    double m_OecX;
    double m_OecY;
    double m_OecZ;
    double m_OecTotCharge;
    double m_OecEnergy;

    TTree* m_ntuple;

    int m_PromptEvt;
    double m_PromptX;
    double m_PromptY;
    double m_PromptZ;
    double m_PromptCharge;

    int m_DelayEvt;
    double m_DelayX;
    double m_DelayY;
    double m_DelayZ;
    double m_DelayCharge;


};


#endif
