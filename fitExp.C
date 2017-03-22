
void fit(double min, double max, TString name1, TString name2="", TString name3="") {
	TChain* ch = new TChain("tree");
	ch->Add(name1.Data());
	if(name2!="") {
		ch->Add(name2.Data());
	}
	if(name3!="") {
		ch->Add(name3.Data());
	}
	ch->GetEntry(0);
	cout << ch->GetLeaf("T0")->GetValue() << endl;
	TString toplot("RateBoard5/1e3 : (T0-");
	toplot += ch->GetLeaf("T0")->GetValue();
	toplot += "+TimeStamp*1/64e6)/60.";
	cout << toplot << endl;
	ch->Draw(toplot.Data(), "Evt != 0 && Evt%50==0");
	TGraph *g = new TGraph(ch->GetSelectedRows(),ch->GetV2(),ch->GetV1());
	g->SetMarkerSize(0.5);
	g->Draw("ap");
	
	TF1* f = new TF1("f", "[0] + [1]*TMath::Exp(-1*[2]*x) + [3]*TMath::Exp(-1*[4]*x)");
	//f->FixParameter(0,0);
	f->SetParameter(1,10000);
	f->FixParameter(2,0.034); // 11C
	//f->FixParameter(4,2.156); // 10C
 	f->FixParameter(4,0.34022); // 15O
	//f->FixParameter(4,0);
	f->SetLineColor(kRed);
g->Fit("f", "", "", min,max);
	cout << "Period1 = " << TMath::Log(2)/f->GetParameter(2) << " minutes" << endl;
	
	cout << "Period2 = " << TMath::Log(2)/f->GetParameter(4) << " minutes" << endl;
}

void fitExp()
{
	/*
	fit(23.6, 140, "~/godaq_rootfiles/analysis_v2.11-calibG2/run83Mult2.root"
 	       ,"~/godaq_rootfiles/analysis_v2.11-calibG2/run84Mult2.root"
// 	       ,"~/godaq_rootfiles/analysis_v2.11-calibG2/run85Mult2.root"
	);*/
	
	
	fit(48.5, 180, "~/godaq_rootfiles/analysis_v2.11-calibG2/run91Mult2.root"
 	       ,"~/godaq_rootfiles/analysis_v2.11-calibG2/run92Mult2.root"
// 	       ,"~/godaq_rootfiles/analysis_v2.11-calibG2/run85Mult2.root"
	);
}