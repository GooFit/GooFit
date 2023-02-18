{
    auto file = new TFile("MC/MC.root");
    tree = (TTree*)file->Get("DecayTree");

    TCanvas foo;
    tree->Draw("thetaprime:mprime","","");
    foo.SaveAs("MC/SqDP.png");
    tree->Draw("mprime","","colz");
    foo.SaveAs("MC/mprime.png");
    tree->Draw("thetaprime","","colz");
    foo.SaveAs("MC/thetaprime.png");
    tree->Draw("s23:s13","","");
    foo.SaveAs("MC/DP.png");
    tree->Draw("s23","","");
    foo.SaveAs("MC/s23.png");
    tree->Draw("s13","","");
    foo.SaveAs("MC/s13.png");
    tree->Draw("s12","","");
    foo.SaveAs("MC/s12.png");
}