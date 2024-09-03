import math

import ROOT
from ROOT import TLorentzVector
from ROOT import TVector3
import uproot

def collect(filename):
    with uproot.open(filename) as f:
        tree = f["events"]
        L_E1_data = tree.arrays(["(LCandidates_h1p**2 + LCandidates_h1m**2)**0.5"])["(LCandidates_h1p**2 + LCandidates_h1m**2)**0.5"]
        L_P1x_data = tree.arrays(["LCandidates_h1px"])["LCandidates_h1px"]
        L_P1y_data = tree.arrays(["LCandidates_h1py"])["LCandidates_h1py"]
        L_P1z_data = tree.arrays(["LCandidates_h1pz"])["LCandidates_h1pz"]
        L_m1_data = tree.arrays(["LCandidates_h1m"])["LCandidates_h1m"]
        L_E2_data = tree.arrays(["(LCandidates_h2p**2 + LCandidates_h2m**2)**0.5"])["(LCandidates_h2p**2 + LCandidates_h2m**2)**0.5"]
        L_P2x_data = tree.arrays(["LCandidates_h2px"])["LCandidates_h2px"]
        L_P2y_data = tree.arrays(["LCandidates_h2py"])["LCandidates_h2py"]
        L_P2z_data = tree.arrays(["LCandidates_h2pz"])["LCandidates_h2pz"]
        L_Q1_data = tree.arrays(["LCandidates_h1q"])["LCandidates_h1q"]
        L_m2_data = tree.arrays(["LCandidates_h2m"])["LCandidates_h2m"]
        Dilepton_mass_data = tree.arrays(["LCandidates_mass**2"])["LCandidates_mass**2"] #tree.arrays(["(LCandidates_h1p**2 + LCandidates_h1m**2) + (LCandidates_h2p**2 + LCandidates_h2m**2) - (LCandidates_h1px + LCandidates_h2px)**2 - (LCandidates_h1py + LCandidates_h2py)**2 - (LCandidates_h1pz + LCandidates_h2pz)**2"])["(LCandidates_h1p**2 + LCandidates_h1m**2) + (LCandidates_h2p**2 + LCandidates_h2m**2) - (LCandidates_h1px + LCandidates_h2px)**2 - (LCandidates_h1py + LCandidates_h2py)**2 - (LCandidates_h1pz + LCandidates_h2pz)**2"]

        Mu_E1_data = tree.arrays(["(Muons_mu1p**2 + Muons_mu1m**2)**0.5"])["(Muons_mu1p**2 + Muons_mu1m**2)**0.5"]
        Mu_P1x_data = tree.arrays(["Muons_mu1px"])["Muons_mu1px"]
        Mu_P1y_data = tree.arrays(["Muons_mu1py"])["Muons_mu1py"]
        Mu_P1z_data = tree.arrays(["Muons_mu1pz"])["Muons_mu1pz"]
        Mu_q1_data = tree.arrays(["Muons_mu1q"])["Muons_mu1q"]

        Mu_E2_data = tree.arrays(["(Muons_mu2p**2 + Muons_mu2m**2)**0.5"])["(Muons_mu2p**2 + Muons_mu2m**2)**0.5"]
        Mu_P2x_data = tree.arrays(["Muons_mu2px"])["Muons_mu2px"]
        Mu_P2y_data = tree.arrays(["Muons_mu2py"])["Muons_mu2py"]
        Mu_P2z_data = tree.arrays(["Muons_mu2pz"])["Muons_mu2pz"]
        Mu_q2_data = tree.arrays(["Muons_mu2q"])["Muons_mu2q"]

    costhetak_array = []
    costhetal_array = []
    phi_array = []
    original = []
    new = []
    for i in range(len(L_E1_data)):
        #if Dilepton_mass_data[i][0] >= 15 and Dilepton_mass_data[i][0] <= 19:
        if Mu_q1_data[i][0] > 0:
            muplus = TLorentzVector(Mu_P1x_data[i][0], Mu_P1y_data[i][0], Mu_P1z_data[i][0], Mu_E1_data[i][0])
            muminus = TLorentzVector(Mu_P2x_data[i][0], Mu_P2y_data[i][0], Mu_P2z_data[i][0], Mu_E2_data[i][0])
        else:
            muminus = TLorentzVector(Mu_P1x_data[i][0], Mu_P1y_data[i][0], Mu_P1z_data[i][0], Mu_E1_data[i][0])
            muplus = TLorentzVector(Mu_P2x_data[i][0], Mu_P2y_data[i][0], Mu_P2z_data[i][0], Mu_E2_data[i][0])
        if L_m1_data[i][0] > 0.5:
            proton = TLorentzVector(L_P1x_data[i][0], L_P1y_data[i][0], L_P1z_data[i][0], L_E1_data[i][0])
            Lzero = L_Q1_data[i][0] > 0
            pion = TLorentzVector(L_P2x_data[i][0], L_P2y_data[i][0], L_P2z_data[i][0], L_E2_data[i][0])
        else:
            pion = TLorentzVector(L_P1x_data[i][0], L_P1y_data[i][0], L_P1z_data[i][0], L_E1_data[i][0])
            Lzero = L_Q1_data[i][0] < 0
            proton = TLorentzVector(L_P2x_data[i][0], L_P2y_data[i][0], L_P2z_data[i][0], L_E2_data[i][0])
        mumu = muminus + muplus
        protonpi = proton + pion
        b = muminus + muplus + proton + pion

        mumuboost = TVector3(-mumu.BoostVector())
        protonpiboost = TVector3(-protonpi.BoostVector())
        bboost = TVector3(-b.BoostVector())

        muminusd = TLorentzVector(muminus)
        muminusd.Boost(mumuboost)
        muplusd = TLorentzVector(muplus)
        muplusd.Boost(mumuboost)
        bd = TLorentzVector(b)
        bd.Boost(mumuboost)
        #costhetal_array.append(math.cos(muplusd.Vect().Angle(-bd.Vect())))
        if (Lzero):
            costhetal_array.append(math.cos(muplusd.Vect().Angle(-bd.Vect())))
        else:
            costhetal_array.append(math.cos(muminusd.Vect().Angle(-bd.Vect())))
        protondd = TLorentzVector(proton)
        protondd.Boost(protonpiboost)
        bdd = TLorentzVector(b)
        bdd.Boost(protonpiboost)
        costhetak_array.append(math.cos(protondd.Vect().Angle(-bdd.Vect())))

        protonddd = TLorentzVector(proton)
        protonddd.Boost(bboost)
        pionddd = TLorentzVector(pion)
        pionddd.Boost(bboost)
        muminusddd = TLorentzVector(muminus)
        muminus.Boost(bboost)
        muplusddd = TLorentzVector(muplus)
        muplusddd.Boost(bboost)

        normalppi = protonddd.Vect().Cross(pionddd.Vect())
        normalmumu = muplusddd.Vect().Cross(muminusddd.Vect())

        protonpiddd = TLorentzVector(protonpi)
        protonpiddd.Boost(bboost)

        if Lzero:
            phi = normalppi.Angle(normalmumu)
            if (normalmumu.Cross(normalppi).Dot(protonpiddd.Vect()) < 0.0):
                phi = -phi
        else:
            phi = normalppi.Angle(-normalmumu)
            if (normalmumu.Cross(normalppi).Dot(protonpiddd.Vect()) < 0.0):
                phi = -phi


        protonddd = TLorentzVector(proton)
        protonddd.Boost(bboost)
        pionddd = TLorentzVector(pion)
        pionddd.Boost(bboost)
        muminusddd = TLorentzVector(muminus)
        muminus.Boost(bboost)
        muplusddd = TLorentzVector(muplus)
        muplusddd.Boost(bboost)
        plminus = TVector3(muminusddd.Vect())
        plplus = TVector3(muplusddd.Vect())
        pproton = TVector3(protonddd.Vect())
        ppion = TVector3(pionddd.Vect())
        #print('sum', (plminus + plplus + ppion + pproton).Mag())
        # print(type(plminus.Cross(plplus)))
        # print(type((1/(plminus.Cross(plplus)).Mag())*(plminus.Cross(plplus))))
        el = (1/(plminus.Cross(plplus)).Mag()) * (plminus.Cross(plplus))
        ek = (1/(pproton.Cross(ppion)).Mag()) * (pproton.Cross(ppion))
        ez = (1/(pproton + ppion).Mag()) * (pproton + ppion)
        # print("mag ", ez.Mag())
        # print(ek.Dot(el)**2 + ((el.Cross(ek)).Dot(ez))**2)
        #print(math.cos(phi), math.sin(phi), ek.Dot(el), ((el.Cross(ek)).Dot(ez)))
        # print("a")
        cos = ek.Dot(el)
        sin = ((el.Cross(ek)).Dot(ez))
        if (sin**2 + cos**2 - 1 > 10e-6):
            print('greater')
        if sin < 0:
            p = math.acos(cos) - math.pi
        else:
            p = math.acos(cos)
        #print("diff ", p - phi)
        #phi_array.append(p)
        # phi_array.append(p)
        phi_array.append(p)
        original.append(phi)
    return (costhetal_array, costhetak_array, phi_array)



















