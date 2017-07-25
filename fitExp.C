using namespace RooFit;

Double_t LambdaToPeriod(Double_t lambda) {
	return TMath::Log(2)/lambda;
}

Double_t exp(Double_t *x,Double_t *par) {
      Double_t fitval = par[0] + par[1]*par[2]/60.*TMath::Exp(-1*par[2]*x[0]);
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

// par[0]: lambda
// par[1]: time spill
// par[2]: rate
// par[3]: constant (for gamma prompt contribution)
Double_t funcBuildUp(Double_t* x, Double_t* par) {
	double ak = 1;
	double RFperiodsec = 1/24.85e6;
	double RFperiodmin = RFperiodsec/60.;
	double k = x[0]/RFperiodmin; // 24.85e6 is the RF frequency
	//cout << "k=" << k << endl;
// 	for(int i = 1; i < k-1; i++) {
// 		cout << "i= " << i << endl;
// 		double temp=1;
// 		//for(int j=0; j<i; j++) {
// 			temp*=TMath::Exp(-1*i*par[0]*(par[1]+par[2]));
// 		//}
// 		ak+=temp;
// 	}
	ak*=(1-TMath::Exp(-1*k*par[0]*RFperiodmin))/(1-TMath::Exp(-1*par[0]*RFperiodmin));
	ak*=(1-TMath::Exp(-1*par[0]*par[1]));
	ak*=TMath::Exp(-1*par[0]*(RFperiodmin - par[1]));
	ak*=par[2];
	ak+=par[3];
	return ak;
}

Double_t expPlusBuildUp(Double_t* x, Double_t* par) {
	Double_t* parExp = new Double_t[3];
	parExp[0] = par[0];
	parExp[1] = par[1];
	parExp[2] = par[2];
	Double_t exp1 = exp(x, parExp);
	
	Double_t* parBuildUp = new Double_t[5];
	parBuildUp[0] = par[3];
	parBuildUp[1] = par[4];
	parBuildUp[2] = par[5];
	parBuildUp[3] = par[6];
	Double_t buildUp = funcBuildUp(x, parBuildUp);
	
	if(x[0] < 0.74) {
		return 0;
	}
	else if(x[0] > 23.5) {
		return exp1;
	} 
	return buildUp;
}

Double_t funcConst(Double_t *x, Double_t *par) {
	if(x[0] < par[0] || x[0] > par[1]) {
		return 0;
	}
	return par[2];
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
	ch->Draw(toplot.Data(), "Evt != 0 && Evt%1000==0" && cut, "goff");
	TGraph *g = new TGraph(ch->GetSelectedRows(),ch->GetV2(),ch->GetV1());
	g->SetMarkerSize(0.8);
	g->SetMarkerStyle(4);
	return g;
}

void fit(TF1* f, TCut cut, double xmin, double xmax, double ymin, double ymax, TString text, TString name1, TString name2="", TString name3="") {
	TGraph* g = getGraph(cut , name1, name2, name3);
	g->GetXaxis()->SetRangeUser(0, xmax+5);
	g->GetYaxis()->SetRangeUser(ymin, ymax);
	g->Draw("ap");
	g->GetXaxis()->SetTitle("Time [minutes]");
	g->GetYaxis()->SetTitle("Rate [MHz]");
	g->GetXaxis()->SetTitleSize(0.05);
	g->GetYaxis()->SetTitleSize(0.05);
	g->GetXaxis()->SetTitleOffset(1.25);
	g->GetYaxis()->SetTitleOffset(1.5);
	g->GetXaxis()->SetLabelSize(0.05);
	g->GetYaxis()->SetLabelSize(0.05);
	gPad->SetGridx();
	gPad->SetGridy();
	PutText(0.45, 0.82, kBlack, "LAPD");
	PutText(0.45, 0.75, kBlack, "Protons 65 MeV, I = 5 nA");
	PutText(0.45, 0.68, kBlack, text.Data());
	g->Fit(f->GetName(), "", "", xmin,xmax);
	cout << "Period1 = " << LambdaToPeriod(f->GetParameter(2)) << " minutes" << endl;
	f->SetLineWidth(3);
	gPad->Update();
	
	//TLegend* leg = new TLegend(0.5,0.36,0.83,0.61);
	TLegend* leg = new TLegend(0.5,0.36,0.83,0.59);
	leg->SetBorderSize(0);
	leg->AddEntry(g, "Observed", "p");
	if(f->GetNpar() == 7) {
		leg->AddEntry(f, "Expected (Total)", "l");
		TF1* fexpPlusBuildUp_1 = new TF1("fexpPlusBuildUp_1", expPlusBuildUp, 0, xmax, 7);
		fexpPlusBuildUp_1->SetParameters(f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3),f->GetParameter(4),f->GetParameter(5), 0);
		fexpPlusBuildUp_1->SetLineColor(kBlack);
		fexpPlusBuildUp_1->SetLineStyle(9);
		fexpPlusBuildUp_1->SetLineWidth(2);
		fexpPlusBuildUp_1->Draw("same");
		leg->AddEntry(fexpPlusBuildUp_1, "Expected (^{11}C)", "l");
		TF1* fConstant = new TF1("fConstant", funcConst, xmin, xmax, 3);
		fConstant->SetParameters(0.74, 23.5, f->GetParameter(6));
		fConstant->SetLineColor(kBlack);
		fConstant->SetLineWidth(2);
		fConstant->SetLineStyle(2);
		fConstant->SetNpx(1e3);
		fConstant->Draw("same");
		leg->AddEntry(fConstant, "Expected (prompt)", "l");
	}
	if(f->GetNpar() == 6) {
		TF1* func2exps_1 = new TF1("func2exps_1", exp, xmin, xmax, 3);
		func2exps_1->SetParameters(f->GetParameter(0),f->GetParameter(1),f->GetParameter(2));
		func2exps_1->SetLineColor(kBlack);
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
	func1exp->SetParNames("Constant_c11", "N_C11","Lambda_C11");
	func1exp->FixParameter(2, 0.034);
	func1exp->SetParLimits(0,0,1000);
	func1exp->SetLineColor(kBlack);
	func1exp->SetLineStyle(3);
	
	xmin = 0;
	TF1* fexpPlusBuildUp = new TF1("fexpPlusBuildUp", expPlusBuildUp, 0, xmax, 7);
	fexpPlusBuildUp->SetNpx(1000);
	fexpPlusBuildUp->SetParLimits(0,0,1000);
	//fexpPlusBuildUp->FixParameter(0,0);
	fexpPlusBuildUp->SetParameter(1,1843);
	fexpPlusBuildUp->FixParameter(2,0.034);
	fexpPlusBuildUp->FixParameter(3,0.034);
	fexpPlusBuildUp->FixParameter(4,3.7e-9/60.);
	fexpPlusBuildUp->SetParameter(5,10);
	fexpPlusBuildUp->SetParameter(6,0.28);
	fexpPlusBuildUp->SetLineColor(kBlack);

	/*
	fexpPlusBuildUp->GetXaxis()->SetTitleSize(0.06);
	fexpPlusBuildUp->GetYaxis()->SetTitleSize(0.06);
	fexpPlusBuildUp->GetXaxis()->SetTitleOffset(1.25);
	fexpPlusBuildUp->GetYaxis()->SetTitleOffset(1.1);
	fexpPlusBuildUp->GetXaxis()->SetLabelSize(0.06);
	fexpPlusBuildUp->GetYaxis()->SetLabelSize(0.06);
	*/

	fit(fexpPlusBuildUp, "", xmin, xmax, 0, 0.85, "HDPE target (5#times5#times5 cm^{3})", "~/godaq_rootfiles/analysis_v3.2-calibG2/run83LOR.root"
 	       ,"~/godaq_rootfiles/analysis_v3.2-calibG2/run84LOR.root");
	c1->SaveAs("ActDesactHDPE.png");
	
	/*
	cout << "***** Target PMMA 5*5 cm ***********" << endl;
	TCanvas* c2 = new TCanvas("c2","c2");
	xmin = 48.5;
	xmax = 160;
	TF1* func2exps = new TF1("f", twoexps, xmin, xmax, 6);
	func2exps->SetParNames("Constant_C11", "N_C11","Lambda_C11", "Constant_O15", "N_O15","Lambda_O15");
	func2exps->FixParameter(2, 0.034);
	func2exps->FixParameter(5, 0.34022);
	func2exps->SetParLimits(0,0,1000);
	func2exps->SetParLimits(3,0,1000);
	func2exps->SetLineColor(kRed);
	
	fit(func2exps, "", xmin, xmax, 0, 1.6, "Target PMMA 5#times5 cm", 
	    "~/godaq_rootfiles/analysis_v2.11-calibG2/run91Mult2.root",
	    "~/godaq_rootfiles/analysis_v2.11-calibG2/run92Mult2.root");
	c2->SaveAs("c2.png");
	
	cout << "***** Target PMMA splitted ***********" << endl; 
	TCanvas* c3 = new TCanvas("c3","c3");
	xmin = 30.2;
	xmax = 140;
	fit(func2exps, "", xmin, xmax, 0, 1.8, "Target PMMA (splitted)", "~/godaq_rootfiles/analysis_v2.11-calibG2/run112Mult2.root");
	c3->SaveAs("c3.png");
	*/
}
