{

gSystem->Load("libRooFit") ;
using namespace RooFit ;

// --- Observable ---
RooRealVar deltaM("deltaM","diff (GeV)", 0.14,0.158) ;//0.14,0.158
  
// --- Parameters ---
RooRealVar sigma("sigma","width of gaussian", 0.0007, 0.0001, 0.003);  //0.0007, 0.0001, 0.003 ;
RooRealVar mean("mean","mean of gaussian",0.1455, 0.14, 0.15) ;

RooRealVar sigma2("sigma2","width of gaussian",0.0006, 0.0002, 0.005) ; //0.0006, 0.0002, 0.005
//RooRealVar mean2("mean2","mean of gaussian",0.146, 0.14, 0.15) ;
   
// --- Build Gaussian PDF ---
RooGaussian sign("sign","signal PDF",deltaM,mean,sigma) ;
RooGaussian sign2("sign2","signal PDF",deltaM,mean,sigma2) ;
RooRealVar fr("fr","width of gaussian",0.999, 0.000, 1.00) ;
RooAddPdf sig("sig", "g1+g2", RooArgSet(sign, sign2),fr);

RooRealVar m0("m0", "m0", 0.1396);
RooRealVar p0("p0","p0", 0.0002, 0., 1);  //0.005,0.,1 ;
RooRealVar p1("p1","p1",10.,0.,50.) ; //4.0,0.,50.
RooRealVar p2("p2","p2",4.0,-100.,100.) ; //-4.0,-100.,200.
RooGenericPdf bkg("bkg","(1-exp(-(deltaM-m0)/p0))*(deltaM/m0)^p1 + p2*(deltaM/m0 - 1)",RooArgSet(m0,deltaM,p0,p1,p2)) ;


RooRealVar nsig("nsig","#signal events",50000,0.,100000) ;
RooRealVar nbkg("nbkg","#background events",5000,0.,40000) ;
RooAddPdf model("model","g+a",RooArgList(sig,bkg),RooArgList(nsig,nbkg)) ;

TFile *f00 = new TFile("MC_4b_arric3.root");
TString name = "input_file";
TH1F *da= (TH1F*) f00->Get(name);

TH1F *da1=(TH1F*)da->Clone();

da1->Sumw2();

da1->Rebin(4);

RooDataHist *data = new RooDataHist("data","diff mass",deltaM,da1);
RooFitResult* r = model.fitTo(*data,Save()) ;
FILE *myfile;
myfile = fopen ("results/input_file.txt", "a");
fprintf(myfile, "N4b_MC %f %f \n", nsig.getValV(), nsig.getError());
fclose(myfile);
r->Print("") ;

RooPlot* mesframe = deltaM.frame() ;
data->plotOn(mesframe) ;
model.plotOn(mesframe) ;
model.plotOn(mesframe,Components(bkg),LineStyle(kDashed)) ;
TCanvas *c = new TCanvas;
TImage *img = TImage::Create();

mesframe->GetXaxis()->SetTitle("M(k#pi#pi#pi#pi_{s}) - M(k#pi#pi#pi) [GeV/c^{2}]");
mesframe->SetTitle("");
mesframe->Draw()  ;

TLatex latexLabel2;
latexLabel2.SetTextSize(0.04);
latexLabel2.SetTextFont(32);
latexLabel2.SetNDC();
latexLabel2.DrawLatex(0.3, 0.93, " CMS Simulation, pp #sqrt{s} = 13 TeV, 2016");

img->FromPad(c);
img->WriteImage("Immagini/"+name+"_4b_MC.png");

gSystem->Exit(0);
}
