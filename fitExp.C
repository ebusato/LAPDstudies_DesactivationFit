using namespace RooFit;

Double_t LambdaToPeriod(Double_t lambda) {
	return TMath::Log(2)/lambda;
}

Double_t exp(Double_t *x,Double_t *par) {
      Double_t fitval = par[0] + par[1]*TMath::Exp(-1*par[2]*x[0]);
      return fitval;
}

Double_t twoexps(Double_t *x, Double_t *par) {
	Double_t* par1 = new Double_t[3];
	par1[0] = par[0];
	par1[1] = par[1];
	par1[2] = par[2];
	Double_t exp1 = exp(x, par1);
	
	Double_t* par2 = new Double_t[3];
	par2[0] = par[3];
	par2[1] = par[4];
	par2[2] = par[5];
	Double_t exp2 = exp(x, par2);
	
	return exp1+exp2;
}

TGraph* getGraph(TCut cut, TString name1, TString name2="", TString name3="") {
	TChain* ch = new TChain("tree");
	ch->Add(name1.Data());
	if(name2!="") {
		ch->Add(name2.Data());
	}
	if(name3!="") {
		ch->Add(name3.Data());
	}
	TString toplot("RateLvsL3/1e6 : (T0-");
	//TString toplot("RateBoard4/1e3 : (T0-");
	ch->GetEntry(0);
	Double_t initTime = ch->GetLeaf("T0")->GetValue();
	toplot += initTime;
	toplot += "+TimeStamp*1/64e6)/60.";
	cout << toplot << endl;
	ch->Draw(toplot.Data(), "Evt != 0 && Evt%1==0" && cut, "goff");
	TGraph *g = new TGraph(ch->GetSelectedRows(),ch->GetV2(),ch->GetV1());
	g->SetMarkerSize(0.8);
	
	return g;
}

void fit(TF1* f, TCut cut, double xmin, double xmax, double ymin, double ymax, TString text, TString name1, TString name2="", TString name3="") {
	TGraph* g = getGraph(cut , name1, name2, name3);
	g->GetXaxis()->SetRangeUser(0, xmax+5);
	g->GetYaxis()->SetRangeUser(ymin, ymax);
	g->Draw("ap");
	g->GetXaxis()->SetTitle("time (minutes)");
	g->GetYaxis()->SetTitle("rate (MHz)");
	gPad->SetGridx();
	gPad->SetGridy();
	PutText(0.5, 0.87, kBlack, "LAPD");
	PutText(0.5, 0.8, kBlack, "Protons 65 MeV, I = 4.5 nA");
	PutText(0.5, 0.73, kBlack, text.Data());
	g->Fit("f", "", "", xmin,xmax);
	cout << "Period1 = " << LambdaToPeriod(f->GetParameter(2)) << " minutes" << endl;
	f->SetLineWidth(3);
	TLegend* leg = new TLegend(0.55,0.45,0.85,0.65);
	leg->SetBorderSize(0);
	leg->AddEntry(g, "Observed", "pl");
	if(f->GetNpar() == 3) {
		leg->AddEntry(f, "Expected (^{11}C)", "l");
		leg->SetY1(0.55);
	}
	if(f->GetNpar() == 6) {
		TF1* func2exps_1 = new TF1("func2exps_1", exp, xmin, xmax, 3);
		func2exps_1->SetParameters(f->GetParameter(0),f->GetParameter(1),f->GetParameter(2));
		func2exps_1->SetLineColor(kBlue);
		func2exps_1->SetLineWidth(2);
		func2exps_1->Draw("same");
		TF1* func2exps_2 = new TF1("func2exps_2", exp, xmin, xmax, 3);
		func2exps_2->SetParameters(f->GetParameter(3),f->GetParameter(4),f->GetParameter(5));
		func2exps_2->SetLineColor(kGreen+2);
		func2exps_2->SetLineWidth(2);
		func2exps_2->Draw("same");
		func2exps_2->Draw("same");
		cout << "Period2 = " << LambdaToPeriod(f->GetParameter(5)) << " minutes" << endl;
		//f->Draw("same");
		leg->AddEntry(f, "Expected (Total)", "l");
		leg->AddEntry(func2exps_1, "Expected (^{11}C)", "l");
		leg->AddEntry(func2exps_2, "Expected (^{15}O)", "l");
	}
	leg->Draw();
	
	return g;
}

void fitExp()
{
	cout << "***** Target HDPE 5*5 cm ***********" << endl;
	TCanvas* c1 = new TCanvas("c1","c1");
	double xmin = 23.6;
	double xmax = 140;
	TF1* func1exp = new TF1("f", exp, xmin, xmax, 3);
	func1exp->SetParNames("Constant", "Norm","Lambda");
	func1exp->FixParameter(2, 0.034);
	func1exp->SetLineColor(kRed);
	fit(func1exp, "", xmin, xmax, 0, 0.85, "Target HDPE 5#times5 cm", "~/godaq_rootfiles/analysis_v2.11-calibG2/run83Mult2.root"
 	       ,"~/godaq_rootfiles/analysis_v2.11-calibG2/run84Mult2.root");
	
	cout << "***** Target PMMA 5*5 cm ***********" << endl;
	TCanvas* c2 = new TCanvas("c2","c2");
	xmin = 48.5;
	xmax = 160;
	TF1* func2exps = new TF1("f", twoexps, xmin, xmax, 6);
	func2exps->SetParNames("Constant1", "Norm1","Lambda1", "Constant2", "Norm2","Lambda2");
	func2exps->FixParameter(2, 0.034);
	func2exps->FixParameter(5, 0.34022);
	func2exps->SetParLimits(0,0,1000);
	func2exps->SetParLimits(3,0,1000);
	func2exps->SetLineColor(kRed);
	
	fit(func2exps, "", xmin, xmax, 0, 1.6, "Target PMMA 5#times5 cm", 
	    "~/godaq_rootfiles/analysis_v2.11-calibG2/run91Mult2.root",
	    "~/godaq_rootfiles/analysis_v2.11-calibG2/run92Mult2.root");
	
	
	cout << "***** Target PMMA splitted ***********" << endl; 
	TCanvas* c3 = new TCanvas("c3","c3");
	xmin = 30.2;
	xmax = 140;
	fit(func2exps, "", xmin, xmax, 0, 1.8, "Target PMMA (splitted)", "~/godaq_rootfiles/analysis_v2.11-calibG2/run112Mult2.root");
}