# -*- coding: utf-8 -*-
"""
Created on Mon Jan  1 13:32:09 2018

@author: hkaveh
"""

import nuSolutions as nu
import ROOT as r
import lhapdf
import numpy as np
import math
from operator import itemgetter
print lhapdf.version()
lhapdf.pathsPrepend("/home/hkaveh/")
print lhapdf.availablePDFSets()
print lhapdf.paths()

pset = lhapdf.getPDFSet('NNPDF23_nlo_as_0118')
print pset.description
p = lhapdf.mkPDF('NNPDF23_nlo_as_0118')
def get_dalitz_prob(lep,top,mb,mw):
    mte = lep.Dot( top );
    mt = top.M();
    mt2 = mt * mt;
    mb2 = mb * mb;
    mw2 = mw * mw;
    mt2_mb2 = mt2 - mb2;

    return 4. * mt*lep.E() * ( mt2 - mb2 - 2. * mt*lep.E() ) /( mt2_mb2 * mt2_mb2 + mw2 * ( mt2 - mb2 ) - 2. * mw2 * mw2 )

def getweight(lep_p,lep_m,nu1,nu2,bquark1,bquark2):
    
    t1 = lep_p + nu1 + bquark1;
    
    t2 = lep_m + nu2 + bquark2;
    e_com = 13000
    top_mass = 172.5
    mw=80.4
    mb=4.8
    
    prob_dalitz = 1.0;
    prob_dalitz *= get_dalitz_prob( lep_p, t1, mb, mw );
    prob_dalitz *= get_dalitz_prob( lep_m, t2, mb, mw );

    #Determine x1 and x2
    x1 = ( t1.E() + t2.E() + t1.Pz() + t2.Pz() ) / e_com;
    x2 = ( t1.E() + t2.E() - t1.Pz() - t2.Pz() ) / e_com;
    
  

    sbar1 = p.xfxQ(-3,x1, top_mass)
    sbar2 = p.xfxQ(-3,x2, top_mass)
    ubar1 = p.xfxQ(-2,x1, top_mass)
    ubar2 = p.xfxQ(-2,x2, top_mass)
    dbar1 = p.xfxQ(-1,x1, top_mass)
    dbar2 = p.xfxQ(-1,x2, top_mass)
    g1    = p.xfxQ(21,x1, top_mass)
    g2    = p.xfxQ(21,x2, top_mass)
    d1    = p.xfxQ(1,x1, top_mass)
    d2    = p.xfxQ(1,x2, top_mass)
    u1    = p.xfxQ(2,x1, top_mass)
    u2    = p.xfxQ(2,x2, top_mass)
    s1    = p.xfxQ(3,x1, top_mass)
    s2    = p.xfxQ(3,x2, top_mass)

    #Should glue-glue be doubled? Probably not, but plot histo later
    pdf_prob = (u1*ubar2 + u2*ubar1 + d1*dbar2 + d2*dbar1 + s1*sbar2 + s2*sbar1 +  g1*g2)

    #print pdf_prob ,' weight ' , prob_dalitz
    return pdf_prob*prob_dalitz

def getKey(item):
    return item[0][0]

f = r.TFile.Open("ElMuDATA.root_tree.root")
f_hist = r.TFile.Open("resolution.root")
tree = r.gDirectory.Get('')
t = f.Get("tree")
rand3 = r.TRandom3()
h_ptResolution = f_hist.Get('ptResolution')
h_phiResolution = f_hist.Get('phiResolution')
print h_ptResolution.GetRandom()
entries = t.GetEntriesFast()
outFile = r.TFile.Open("solutionsTemp.root","recreate")
outTree = r.TTree( 't1' , '')
# create 1 dimensional float arrays (python's float datatype corresponds to c++ doubles)
# as fill variables
lepton1 = np.zeros(1, dtype=float)
lepton2 = np.zeros(1, dtype=float)


top1 = np.zeros(1, dtype=float)
top2 = np.zeros(1, dtype=float)
cos1 = np.zeros(1, dtype=float)
cos2 = np.zeros(1, dtype=float)
mtt = np.zeros(1, dtype=float)
mtt[0] = -999
cos1[0] = -999
cos2[0] = -999
lepton1[0]=-999
lepton2[0]=-999
top1[0]=-999
top2[0]=-999
#create histogram for cos
h_cos = r.TH1D("cos","",100,-1,1)
# create the branches and assign the fill-variables to them
outTree.Branch('lepton1', lepton1, 'lepton1/D')
outTree.Branch('top1', top1, 'top1/D')
outTree.Branch('lepton2', lepton1, 'lepton2/D')
outTree.Branch('top2', top2, 'top2/D')
outTree.Branch('cos1',cos1,'cos1/D')
outTree.Branch('cos2',cos2,'cos2/D')
outTree.Branch('mtt',mtt,'mtt/D')
for jentry in xrange(20000):
    if jentry % 1000 == 0:
        print "Event", jentry, "/", entries
    t.GetEntry(jentry)
    bList = []
    lList = []
    metxList = []
    metyList = []
    sumjetpx = 0.
    sumjetpy = 0.
    sumleppx = 0.
    sumleppy = 0.
    for pt,eta,phi,e,ptR,phiR,sf in map(None,t.pt_bJets,t.eta_bJets,t.phi_bJets,t.e_bJets,t.ptRes_bJets,t.phiRes_bJets,t.sf_bJets):
        b = r.TLorentzVector()
        b.SetPtEtaPhiE(pt,eta,phi,e)
        bb = b,ptR,phiR,sf
        bList.append(bb)
        sumjetpx+=b.Px()
        sumjetpy+=b.Py()
        
    for pt,eta,phi,e,charge in map(None,t.pt_Leptons,t.eta_Leptons,t.phi_Leptons,t.e_Leptons,t.charge_Leptons):
        l = r.TLorentzVector()
        l.SetPtEtaPhiE(pt,eta,phi,e)
        ll= l,charge
        lList.append(ll)
        sumleppy+=l.Px()
        sumleppx+=l.Py()
    for px,py in map(None,t.px_mets,t.py_mets):
        metxList.append(px)
        metyList.append(py)
    t_weight = []
    for i in xrange(2):
        print i
        if(len(bList) >= 2):
            if i == 0 :
                b1 = bList[0]
                b2 = bList[1]
                #print b1[0]
            elif i == 1:
                b1 = bList[1]
                b2 = bList[0]
                #print b1
            
            
            if(lList[0][1] > 0 ):
                l1 = lList[0][0]
                l2 = lList[1][0]
            else:
                l1 = lList[1][0]
                l2 = lList[0][0]
                
            c = r.TMath.Cos(b1[0].Angle(l1.Vect()))
            print "Event", jentry, "/", entries
                
                
            for smear in xrange(10):
                b1smear = r.TLorentzVector()
                b2smear = r.TLorentzVector()
                #ptNormaldist = r.TRandom.Gaus(0,h_ptResolution.GetRandom())
                #phiNormaldist = r.TRandom.Gaus(0,h_phiResolution.GetRandom())
                b1sf = b1[3]
                b1ptResG = rand3.Gaus(0,b1[1]*math.sqrt(max(b1sf*b1sf-1,0)))
                b1phiResG = rand3.Gaus(0,b1[2])
                
                b1smfac = 1 + b1ptResG
                #print b1smfac,b1sf,b1ptResG
                b2sf = b2[3]
                b2ptResG = rand3.Gaus(0,b2[1]*math.sqrt(max(b1sf*b1sf-1,0)))
                b2phiResG = rand3.Gaus(0,b2[2])
                
                b2smfac = 1 + b2ptResG
                    
                    
                b1smear.SetPtEtaPhiE(b1[0].Pt()*b1smfac,b1[0].Eta(),b1[0].Phi(),b1[0].E())
                b2smear.SetPtEtaPhiE(b2[0].Pt()*b2smfac,b2[0].Eta(),b2[0].Phi(),b2[0].E())
                bsm1 = b1smear
                bsm2 = b2smear
                    
                    
                metx = metxList[0]
                mety = metyList[0]
                unclust_metx = metx + sumjetpx + sumleppx
                unclust_mety = mety + sumjetpy + sumleppy
                metG = rand3.Gaus(0,0.1)
                un_metx_sm = unclust_metx*(1+metG)
                un_mety_sm = unclust_mety*(1+metG)
                metx_sm = metx + sumjetpx - unclust_metx - b1smear.Px()-b2smear.Px()+un_metx_sm
                mety_sm = mety + sumjetpy - unclust_mety - b1smear.Py()-b2smear.Py()+un_mety_sm
                    
                #print b1smear.Pt(),'  ', b2smear.Pt()
                #print b1smear.Phi(),' ', b2smear.Phi()
                dum_neu = r.TLorentzVector()
                try:
                    solver = nu.doubleNeutrinoSolutions((bsm1,bsm2),(l1,l2),(metx_sm,mety_sm),(nu.mW)**2,(nu.mT)**2)
                    
                    Neu = solver.nunu_s
                    Ne = tuple
                    Neutrinos = []
                    for Ne in Neu:
                        neu = r.TLorentzVector()
                        neu.SetXYZM(Ne[0][0],Ne[0][1],Ne[0][2],0)
                        Neutrinos.append(neu)
                        #print Ne[0][2]
                        if(len(Neutrinos) == 2):
                            weight = getweight(l1,l2,Neutrinos[0],Neutrinos[1],bsm1,bsm2)
                            dum_t = weight,Neutrinos[0],Neutrinos[1]
                            t_weight.append(dum_t)
                                
                                
                                
                        if(len(Neutrinos) == 4):
                            weight1 = getweight(l1,l2,Neutrinos[0],Neutrinos[1],bsm1,bsm2)
                            weight2 = getweight(l1,l2,Neutrinos[2],Neutrinos[3],bsm1,bsm2)
                            dum_t1 = weight1,Neutrinos[0],Neutrinos[1]
                            t_weight.append(dum_t1)
                            dum_t2 = weight2,Neutrinos[2],Neutrinos[3]
                            t_weight.append(dum_t2)
                        
                    
                                
                                
                except :
                    continue 
                

    dummy = 0
    highest_item = 0
    if(len(t_weight)>0):
        for item in xrange(len(t_weight)):
            if(math.fabs(t_weight[item][0]) > dummy):
                dummy = math.fabs(t_weight[item][0])
                highest_item = item
                result = t_weight[highest_item]
        if (len(result)>0):
            print result[0]
            if(lList[0][1] > 0 ):
                l1 = lList[0][0]
                l2 = lList[1][0]
            else:
                l1 = lList[1][0]
                l2 = lList[0][0]
            t1 = r.TLorentzVector()
            t2 = r.TLorentzVector()
            w1 = r.TLorentzVector()
            w2 = r.TLorentzVector()
            w1 = (l1 + result[1])
            w2 = (l2 + result[2])
            t1 = (l1+result[1]+bList[0][0])
            t2 = (l2+result[2]+bList[1][0])
            
            tt = t1+t2
            
            ww1 = w1
            ww2 = w2
            ll1 = l1
            ll2 = l2
            
            lepton1[0] = l1.Pt()
            top1[0] = t1.Pt()
            lepton2[0] = l2.Pt()
            top2[0] = t2.Pt()
            mtt[0] = tt.M()
           
            ll1.Boost(-w1.BoostVector())
            ww1.Boost(-t1.BoostVector())
            theta1 = ww1.Angle(ll1.Vect())
            cos1[0] = r.TMath.Cos(theta1)
            h_cos.Fill(r.TMath.Cos(theta1))
            ll2.Boost(-w2.BoostVector())
            ww2.Boost(-t2.BoostVector())
            theta2 = ww2.Angle(ll2.Vect())
            cos2[0] = r.TMath.Cos(theta2)
            h_cos.Fill(r.TMath.Cos(theta2))
            outTree.Fill()    
               
                
            
        
    
outFile.cd()
h_cos.Write()
outFile.Write()
outFile.Close()
