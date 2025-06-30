#include "IBDSelectionAlg.h"
#include "TOF.h"
#include "EvtNavigator/NavBuffer.h"
#include "EvtNavigator/EvtNavHelper.h"
#include "SniperKernel/AlgFactory.h"
#include "SniperKernel/SniperDataPtr.h"
#include "SniperKernel/SniperPtr.h"

#include "Event/CdLpmtCalibHeader.h"
#include "Event/CdSpmtCalibHeader.h"
#include "Event/WpCalibHeader.h"

#include "Event/CdTriggerHeader.h" 
#include "Event/CdTriggerEvt.h"

#include "Event/CdVertexRecHeader.h"

#include "Event/OecHeader.h"
#include "Event/OecEvt.h"
#include "OECTagSvc/OECTagSvc.h"
#include "OECTagID/OECTagID.h"

#include "Identifier/Identifier.h"
#include "Identifier/CdID.h"
#include "Identifier/WpID.h"
#include "RootWriter/RootWriter.h"
#include "TTree.h"
#include "BufferMemMgr/IDataMemMgr.h"
#include "Geometry/IPMTParamSvc.h"

DECLARE_ALGORITHM(IBDSelectionAlg);

IBDSelectionAlg::IBDSelectionAlg(const std::string& name) : AlgBase(name)
{
    m_iEvt = -1;
	m_DelayIBD = -1;
	m_DelayEvt = -1;
	foundDelay = false;
}

bool IBDSelectionAlg::initialize(){
    //----------------------------------------------------------------------------
	const std::string compilation_date = __DATE__;
	const std::string compilation_time = __TIME__;
	std::cout <<"##################################################################"<<std::endl
	<<"The source file was compiled on " << compilation_date<< " at " << compilation_time <<std::endl
	<<"##################################################################"<<std::endl;
	//----------------------------------------------------------------------------

    // =======================================================================
    // Loading PMT positions
    // =======================================================================
	
	SniperPtr<IPMTParamSvc> pmtsvc(getParent(), "PMTParamSvc");
	if (ALL_LPMT_pos.size()==0 && pmtsvc.valid()) {
		TotalLPMT = pmtsvc->get_NTotal_CD_LPMT();
		
		std::cout << " PMT Information " << std::endl;
		
		for (unsigned int ith = 0; ith < TotalLPMT; ith++)
		{
			TVector3 all_pmtCenter(pmtsvc->getPMTX(ith), pmtsvc->getPMTY(ith), pmtsvc->getPMTZ(ith));
			ALL_LPMT_pos.push_back(all_pmtCenter);
		}
	}
	if (ALL_SPMT_pos.size()==0 && pmtsvc.valid()) {
		TotalSPMT = pmtsvc->get_NTotal_CD_SPMT();
		
		for (unsigned int ith = 0; ith < TotalSPMT; ith++)
		{
			TVector3 all_pmtCenter(pmtsvc->getPMTX(ith+20000), pmtsvc->getPMTY(ith+20000), pmtsvc->getPMTZ(ith+20000));
			ALL_SPMT_pos.push_back(all_pmtCenter);
		}
	}

	// =======================================================================
    // GET EVENT
    // =======================================================================
	
	
	LogDebug << "initializing" << std::endl;
	std::cout<<"36"<<std::endl;
	gDirectory->pwd();
	
	SniperDataPtr<JM::NavBuffer> navBuf(getRoot(),"/Event");
	if ( navBuf.invalid() ) {
		LogError << "cannot get the NavBuffer @ /Event" << std::endl;
		return false;
	}
	m_buf = navBuf.data();
	// std::cout << "Buffer size: " << navBuf->size() <<std::endl;
	
    SniperPtr<OECTagSvc> tagsvc(getParent(), "OECTagSvc");
	if(tagsvc.invalid()){
		LogError << "Unable to locate tagsvc" << std::endl;
		return false;
	}
	m_tagsvc = tagsvc.data();
    i_pBiPo214 = m_tagsvc->getpTag("BiPo214Pair");
	i_dBiPo214 = m_tagsvc->getdTag("BiPo214Pair");


	// =======================================================================
    // GET INTERFACE LEVEL
    // =======================================================================
	
	TFile* f = TFile::Open("/sps/juno/mlecocq/Commissioning/massVolumeLS.root", "READ");
	if (!f || f->IsZombie()) {
		std::cerr << "Error opening file!" << std::endl;
		return false;
    }

	gInterfaceLevel = (TGraph*)f->Get("inter_level_s1_graph");
	if (!gInterfaceLevel) {
        std::cerr << "TGraph not found in file!" << std::endl;
        f->Close();
        return false;
    }

	f->Close();

	bookTree();
	return true;
}

bool IBDSelectionAlg::execute(){

    std::cout << "executing: " << ++m_iEvt << std::endl;
	gDirectory->pwd();

    JM::EvtNavigator* nav = 0;
	JM::CdTriggerEvt *triggerevent = 0;
    JM::CdLpmtCalibEvt* calibeventLPMT = 0;
	JM::WpCalibEvt* wpcalibevt = 0;
    JM::OecEvt* oecevt = 0;

    nav = m_buf->curEvt();

	std::cout << "m_iEvt: " << m_iEvt << std::endl;
	std::cout << "Delay Entry: " << m_DelayIBD << std::endl;

	std::cout<<"71"<<std::endl;
	gDirectory->pwd();
	const auto& paths = nav->getPath();
	const auto& refs = nav->getRef();

	LogInfo << "Detector type is  " << nav->getDetectorType()<<std::endl;
	LogInfo << "Start to Explore SmartRef: " << std::endl;
	LogInfo << "Size of paths: " << paths.size() << std::endl;
	LogInfo << "Size of refs: " << refs.size() << std::endl;

    for (size_t i = 0; i < paths.size(); ++i) {
        LogInfo << refs[i]<<" -> ref: " << std::endl;
        const std::string& path = paths[i];
        JM::SmartRef* ref = refs[i];
        JM::EventObject* evtobj = ref->GetObject();

        LogInfo << " path: " << path
            << " ref->entry(): " << ref->entry()
            << " evtobj: " << evtobj 
			<< std::endl;
	}

    // ===================================================================
	// TRIGGEREVT
	// ===================================================================

	auto triggerheader = JM::getHeaderObject<JM::CdTriggerHeader>(nav);
	if(triggerheader){
		// triggerevent = dynamic_cast<JM::CdTriggerEvt*>(triggerheader->event());
		triggerevent = triggerheader->event();
		std::cout<<" CD TriggerEvent Read in: " << triggerevent <<std::endl;
	}

    // ===================================================================
	// CALIBEVT
	// ===================================================================

	auto calibheaderLPMT = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav);
	if (calibheaderLPMT) {
		// calibeventLPMT = dynamic_cast<JM::CdLpmtCalibEvt*>(calibheaderLPMT->event());
		calibeventLPMT = calibheaderLPMT->event();
		LogInfo << "CalibEventLPMT Read in: " << calibeventLPMT << std::endl;
	}
	auto wpcalibhdr = JM::getHeaderObject<JM::WpCalibHeader>(nav);
    if (wpcalibhdr) {
        wpcalibevt = dynamic_cast<JM::WpCalibEvt*>(wpcalibhdr->event());
		LogInfo << "WpEvent Read in: " << wpcalibevt << std::endl;
    }

    // ===================================================================
	// OECEVT
	// ===================================================================

	auto oecheader = JM::getHeaderObject<JM::OecHeader>(nav);
	if(oecheader){
		oecevt = dynamic_cast<JM::OecEvt*>(oecheader->event("JM::OecEvt"));
		LogInfo << "OecEvent Read in: " << oecevt << std::endl;
	}

	m_TriggerType.clear();
	m_PmtId.clear();
	m_HitTime.clear();

    
    const auto& timestamp = nav->TimeStamp();
    m_iRun = nav->RunID();
    m_TimeStamp = timestamp.GetSec()*1000000000ULL + timestamp.GetNanoSec();

    if(triggerevent){
        const auto& type = triggerevent->triggerType();
        const auto& pmtFired = triggerevent->nHitMultiplicity();
		// const auto& volID = triggerevent->volumeId();
		const auto& trigTime = triggerevent->triggerTime();

		// m_TriggerTime = trigTime.GetSec()*1000000000ULL + trigTime.GetNanoSec();

		std::cout<< "Trigger type size " <<type.size()<<std::endl;
		std::cout<< "Triggered PMT size " <<pmtFired<<std::endl;
		// std::cout<< "volume size " <<volID.size()<<std::endl;
		std::cout<< "Trigger Time size " <<trigTime.GetNanoSec()<<std::endl;

		for(auto it = 0; it<type.size(); it++){
			std::cout<<"Trigger type = "<<type[it]<<std::endl;
			m_TriggerType.push_back(type[it]);

		}
    }

	if(calibeventLPMT && oecevt){

		m_OecTotCharge = oecevt->getTotalCharge();
		std::cout << "OecCharge " << m_OecTotCharge << std::endl;
		m_OecX = oecevt->getVertexX();
		m_OecY = oecevt->getVertexY();
		m_OecZ = oecevt->getVertexZ();
		m_OecEnergy = oecevt->getEnergy();

		const auto& chhlistLPMT = calibheaderLPMT->event()->calibPMTCol();
				// std::cout<<"NLPMT Hitted "<<chhlistLPMT.size()<<std::endl;
		for (auto chit = chhlistLPMT.begin(); chit!=chhlistLPMT.end(); ++chit) {
			auto calib = *chit;

			unsigned int pmtId = calib->pmtId();
			Identifier id = Identifier(pmtId);
			int TruePM=CdID::module(id);
			// std::cout << " - PMTID: " << TruePM << std::endl;

			for(unsigned int j=0;j<calib->size();j++){
				m_HitTime.push_back(calib->time(j));
				m_PmtId.push_back(TruePM);
			}
		}

		if(m_iEvt == m_DelayEvt){
			m_DelayHitTimeTOF.clear();
			for (size_t i = 0; i < m_HitTime.size(); i++){
				double timeTOF = 0.0;
				int pmtid = m_PmtId.at(i);
				TVector3 vertex(m_DelayX, m_DelayY, m_DelayZ);
				double interface = (gInterfaceLevel->Eval(m_TimeStamp * 1e-9) - 17.7)*1e3;

				TOFCalculator TOF(vertex, ALL_LPMT_pos.at(pmtid), interface);

				timeTOF = m_HitTime.at(i) - TOF.CalLTOF();
				m_DelayHitTimeTOF.push_back(timeTOF);
			}
			m_ntuple->Fill();
		}
	}



	bool EventCheck = (m_iEvt - (m_DelayIBD + 1)) == 0;
	std::cout << "Event Check condition " << EventCheck << std::endl;
	// std::cout << "Trigger Type " << m_TriggerType[0] << std::endl;

    if(oecevt && calibeventLPMT && EventCheck && m_TriggerType[0] == "nHit"){
        
        // Select Prompt
        std::cout << "Selecting Prompt " << std::endl;
        
        double r_pos = sqrt( pow(m_OecX,2) + pow(m_OecY,2) + pow(m_OecZ,2) );
        bool r_cut = r_pos < 16500 && m_OecZ < 15500;
        
        if( (m_OecTotCharge >= 3000 && m_OecTotCharge <= 20000 && m_OecEnergy != -1) && r_cut ){
            
            std::cout << "Looking for IBD Signal" << std::endl;
            
            std::vector<JM::OecEvt*> prompt_Evt;
            prompt_Evt.push_back(oecevt);
            
            // Check if Delay
            if(findCorrelation(m_buf->current(), prompt_Evt)){
                std::cout << "Found IBD with dt = " << m_TimeDifference << "X: " << m_DelayX << "Y: " << m_DelayY << "Z: " << m_DelayZ << std::endl;

                m_PromptEvt = m_iEvt;
                m_PromptCharge = oecevt->getTotalCharge();
                m_PromptX = oecevt->getVertexX();
                m_PromptY = oecevt->getVertexY();
                m_PromptZ = oecevt->getVertexZ();

				m_EventTag.clear();
                m_EventTag.push_back("IBD");
                m_DelayType = "IBD";

				m_PromptHitTimeTOF.clear();
                for (size_t i = 0; i < m_HitTime.size(); i++){
                    double timeTOF = 0.0;
                    int pmtid = m_PmtId.at(i);
					TVector3 vertex(m_PromptX, m_PromptY, m_PromptZ);
					double interface = (gInterfaceLevel->Eval(m_TimeStamp * 1e-9) - 17.7) * 1e3;
					TOFCalculator TOF(vertex, ALL_LPMT_pos.at(pmtid), interface);

                    if(pmtid <= 18000){
                        timeTOF = m_HitTime.at(i) - TOF.CalLTOF();
						m_PromptHitTimeTOF.push_back(timeTOF);
                    }
                }

            }
            else{
                std::cout << "No Delay Signal Found " <<std::endl;
            }
        }

	}
    if(!foundDelay){
        ++m_DelayIBD;
    } 

	return true;

}

bool IBDSelectionAlg::finalize(){
	LogDebug << "finalizing" << std::endl;
	return true;
}






bool IBDSelectionAlg::isMuonVetoed(JM::NavBuffer::Iterator navit){

	JM::OecHeader* tHeaderOEC = JM::getHeaderObject<JM::OecHeader>(navit->get());
	JM::OecEvt* tEventOEC = dynamic_cast<JM::OecEvt*>(tHeaderOEC->event("JM::OecEvt"));

	const TTimeStamp& ttime = tEventOEC->getTime();
	
	uint64_t dtime = 0.0;
	bool isMuon = false;
	
	std::cout << "\n Muon Search..." << std::endl;
	for(JM::NavBuffer::Iterator tmpit = --navit; tmpit != m_buf->begin(); --tmpit){

		JM::OecHeader* bHeaderOEC = JM::getHeaderObject<JM::OecHeader>(tmpit->get());
		JM::OecEvt* bEventOEC = dynamic_cast<JM::OecEvt*>(bHeaderOEC->event("JM::OecEvt"));

		const TTimeStamp& beforetime = bEventOEC->getTime();
		dtime = ((ttime.GetSec() -  beforetime.GetSec())*1000000000ULL + (ttime.GetNanoSec() - beforetime.GetNanoSec())) * 1e-3; // in us

		if(tmpit == m_buf->current()){
			continue;
		}
		if(dtime > 2000){
			std::cout << "-- No Muon found!" << std::endl;
			break;
		}
		if(bEventOEC->getTotalCharge() > 50000 && dtime < 2000){
			isMuon = true;
			break;
		}
	}

	std::cout << " No Muon " << !isMuon << std::endl;
	return !isMuon;

}

bool IBDSelectionAlg::isIsolated(JM::NavBuffer::Iterator navit){

	
	JM::OecHeader* tHeaderOEC = JM::getHeaderObject<JM::OecHeader>(navit->get());
	JM::OecEvt* tEventOEC = dynamic_cast<JM::OecEvt*>(tHeaderOEC->event("JM::OecEvt"));
	
	const TTimeStamp& ttime = tEventOEC->getTime();
	uint64_t dtime = 0.0;
	double distance = 0.0;
	
	bool isMultiplicity = false;
	
	std::cout << "\n Multiplicity Search..." << std::endl;
	
	//---------------------- Check Multiplicity 2ms before Delay ---------------------
	for(JM::NavBuffer::Iterator preIt = (navit - 1); preIt != m_buf->begin(); --preIt){

		JM::OecHeader* bHeaderOEC = JM::getHeaderObject<JM::OecHeader>(preIt->get());
		JM::OecEvt* bEventOEC = dynamic_cast<JM::OecEvt*>(bHeaderOEC->event("JM::OecEvt"));

		const TTimeStamp& beforetime = bEventOEC->getTime();
		dtime = ((ttime.GetSec() -  beforetime.GetSec())*1000000000ULL + (ttime.GetNanoSec() - beforetime.GetNanoSec())) * 1e-3; // in us

		distance = sqrt(
			pow(tEventOEC->getVertexX() - bEventOEC->getVertexX(), 2) +
			pow(tEventOEC->getVertexY() - bEventOEC->getVertexY(), 2) +
			pow(tEventOEC->getVertexZ() - bEventOEC->getVertexZ(), 2)
		);

		bool charge_cut = (bEventOEC->getTotalCharge() >= 3000 && bEventOEC->getTotalCharge() <= 20000) && bEventOEC->getEnergy() != 1;

		if(preIt == m_buf->current()){
			continue;
		}
		if(dtime > 2000){
			std::cout << "-- No Multiplicity before !" << std::endl;
			break;
		}
		if(distance < 1500 && dtime < 2000 && charge_cut){
			isMultiplicity = true;
			break;
		}
	}

	
	//---------------------- Check Multiplicity 1ms after Delay ---------------------

	for(JM::NavBuffer::Iterator postIt = (navit + 1); postIt != m_buf->end(); ++postIt){

		std::cout << "Post Delay Nav " << postIt->get() << std::endl; 

		JM::OecHeader* aHeaderOEC = JM::getHeaderObject<JM::OecHeader>(postIt->get());
		JM::OecEvt* aEventOEC = dynamic_cast<JM::OecEvt*>(aHeaderOEC->event("JM::OecEvt"));

		const TTimeStamp& aftertime = aEventOEC->getTime();
		dtime = ((aftertime.GetSec() - ttime.GetSec())*1000000000ULL + (aftertime.GetNanoSec() -  ttime.GetNanoSec())) * 1e-3; // in us

		std::cout << "Time difference " << dtime << std::endl;

		distance = sqrt(
			pow(aEventOEC->getVertexX() - tEventOEC->getVertexX(), 2) +
			pow(aEventOEC->getVertexY() - tEventOEC->getVertexY(), 2) +
			pow(aEventOEC->getVertexZ() - tEventOEC->getVertexZ(), 2)
		);

		std::cout << "Distance " << distance << std::endl;


		bool charge_cut = (aEventOEC->getTotalCharge() > 3700 && aEventOEC->getTotalCharge() < 6000) && aEventOEC->getEnergy() != 1;
		if(dtime > 1000){
			std::cout << "-- No Multiplicity after !" << std::endl;
			break;
		}
		if(distance < 1500 && dtime < 1000 && charge_cut){
			isMultiplicity = true;
			break;
		}
	}

	std::cout << " No Multiplicity " << !(isMultiplicity) << std::endl;
	return !(isMultiplicity);
}


bool IBDSelectionAlg::findCorrelation(JM::NavBuffer::Iterator navit, std::vector<JM::OecEvt*>& pEvt){

	std::cout << "Running IBD selection " << std::endl;
	foundDelay = false;

	const TTimeStamp& ttime = pEvt[0]->getTime();

	int DelayCount = 0;
	int offset = 0;

	unsigned long dtime = 0.0;

	for(JM::NavBuffer::Iterator tmpit = ++navit; tmpit != m_buf->end(); ++tmpit){

		std::cout << "Searching..." << std::endl;
		offset++;

		JM::OecHeader* aHeaderOEC = JM::getHeaderObject<JM::OecHeader>(tmpit->get());
		JM::OecEvt* aEventOEC = dynamic_cast<JM::OecEvt*>(aHeaderOEC->event("JM::OecEvt"));

		const TTimeStamp& afterTime = aEventOEC->getTime();
		double aEnergy = aEventOEC->getEnergy();
		double distance = sqrt(
			pow(aEventOEC->getVertexX() - pEvt.at(0)->getVertexX(), 2) +
			pow(aEventOEC->getVertexY() - pEvt.at(0)->getVertexY(), 2) +
			pow(aEventOEC->getVertexZ() - pEvt.at(0)->getVertexZ(), 2)
		);

		dtime = ((afterTime.GetSec() - ttime.GetSec())*1000000000ULL + (afterTime.GetNanoSec() - ttime.GetNanoSec())) * 1e-3; //in us
		std::cout << "Current event dt from prompt: " << dtime << " us" << std::endl;

		double afterTotalCharge = aEventOEC->getTotalCharge();

		bool dt_cut = dtime > 10 && dtime < 1000;
		bool position_cut = distance < 1500;

		if(dtime > 1000){
			break;
		}
		if(position_cut && dt_cut && (afterTotalCharge > 3700 && afterTotalCharge < 6000) && aEnergy != -1){
			m_DelayCharge = afterTotalCharge;
			m_TimeDifference = dtime * 1e3;
			m_DelayX = aEventOEC->getVertexX();
			m_DelayY = aEventOEC->getVertexY();
			m_DelayZ = aEventOEC->getVertexZ();
			DelayIt = tmpit;
			std::cout << "Delay Nav " << tmpit->get() << std::endl; 
			// bool MuonVeto = isMuonVetoed(DelayIt);
			// foundDelay = isIsolated(DelayIt) && MuonVeto;
			foundDelay = isMuonVetoed(DelayIt) && isIsolated(DelayIt);
			break;
		}
	}

	if(foundDelay){
		m_DelayIBD = m_iEvt + offset;
		m_DelayEvt = m_DelayIBD;
	}
	std::cout << " Offset " << offset << "\n" << std::endl;
	
	return foundDelay;
}

bool IBDSelectionAlg::bookTree(){

    SniperPtr<RootWriter> svc(*getRoot(), "RootWriter");

    m_ntuple = svc->bookTree(*m_par, "Data/Prompt", "Prompt Information");
    m_ntuple->Branch("PromptEntry", &m_PromptEvt, "EntryNb/I");
    m_ntuple->Branch("RunId", &m_iRun, "RunId/I");
    m_ntuple->Branch("Tag", &m_EventTag);
    m_ntuple->Branch("TimeStamp", &m_TimeStamp);
    m_ntuple->Branch("PromptCharge", &m_PromptCharge);
    m_ntuple->Branch("PromptX", &m_PromptX);
    m_ntuple->Branch("PromptY", &m_PromptY);
    m_ntuple->Branch("PromptZ", &m_PromptZ);
    m_ntuple->Branch("PromptHitTimeTOF", &m_PromptHitTimeTOF);
    m_ntuple->Branch("DelayEntry", &m_DelayEvt);
    m_ntuple->Branch("DelayCharge", &m_DelayCharge);
    m_ntuple->Branch("DelayX", &m_DelayX);
    m_ntuple->Branch("DelayY", &m_DelayY);
    m_ntuple->Branch("DelayZ", &m_DelayZ);
    m_ntuple->Branch("DelayHitTimeTOF", &m_DelayHitTimeTOF);
    m_ntuple->Branch("TimeDifference", m_TimeDifference);

	return true;
}
