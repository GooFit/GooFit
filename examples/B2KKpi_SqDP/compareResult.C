{
    auto file = new TFile("MC/MC.root");
    auto tree = (TTree*)file->Get("DecayTree");

    auto file_lau = new TFile("MC/toyKKpi.root");
    auto tree_lau = (TTree*)file_lau->Get("genResults");

    TCanvas foo;
    tree->Draw("thetaprime:mprime","","");
    foo.SaveAs("MC/SqDP.png");
    tree_lau->Draw("mPrime","charge>0.","hist");
    tree->Draw("mprime","","Esame");
    foo.SaveAs("MC/mprime.png");
    tree_lau->Draw("thPrime","charge>0.","hist");
    tree->Draw("thetaprime","","Esame");
    foo.SaveAs("MC/thetaprime.png");
    tree->Draw("s23:s13","","");
    foo.SaveAs("MC/DP.png");

    tree_lau->Draw("m13Sq","charge>0.","hist");
    tree->Draw("s13","","Esame");
    foo.SaveAs("MC/s13.png");

    tree_lau->Draw("m12Sq","charge>0.","hist");
    tree->Draw("s12","","Esame");
    foo.SaveAs("MC/s12.png");

    tree_lau->Draw("m23Sq","charge>0.","hist");
    tree->Draw("s23","","Esame");
    foo.SaveAs("MC/s23.png");
}