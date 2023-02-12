{
    auto file = new TFile("MC/MC.root");
    auto tree = (TTree*)file->Get("DecayTree");

    auto file_lau = new TFile("MC/gen-3K.root");
    auto tree_lau = (TTree*)file_lau->Get("genResults");

    TCanvas foo;
    tree->Draw("thetaprime:mprime","","");
    foo.SaveAs("MC/SqDP.png");
    tree_lau->Draw("mPrime","","hist");
    tree->Draw("mprime","","Esame");
    foo.SaveAs("MC/mprime.png");
    tree_lau->Draw("thPrime","","hist");
    tree->Draw("thetaprime","","Esame");
    foo.SaveAs("MC/thetaprime.png");
    tree->Draw("s23:s13","","");
    foo.SaveAs("MC/DP.png");
    tree_lau->Draw("m13Sq","","hist");
    tree->Draw("s13","","Esame");
    gPad->SetLogy();
    foo.SaveAs("MC/s13.png");
    gPad->SetLogy(0);
    tree_lau->Draw("m12Sq","","hist");
    tree->Draw("s12","","Esame");
    foo.SaveAs("MC/s12.png");
    tree_lau->Draw("m23Sq","","hist");
    tree->Draw("s23","","Esame");
    foo.SaveAs("MC/s23.png");
}