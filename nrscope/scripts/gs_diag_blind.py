import csv,sys,numpy as np,collections as C
g=[r for r in csv.DictReader(open(sys.argv[1]+"/golden_sniffer_grants.csv")) if int(r["blind_pci"])>=0]
c=np.array([float(r["blind_pci_coh"]) for r in g])
print("invisible grants with energy (> 6 dB) searched over all 1008 PCIs: %d"%len(g))
print("best-PCI coherence: p10 %.2f median %.2f p90 %.2f"%tuple(np.percentile(c,[10,50,90])))
good=[r for r in g if float(r["blind_pci_coh"])>0.6]
print("found the DM-RS of another PCI (> 0.6): %d (%.0f%%), PCIs: %s"%(len(good),100*len(good)/max(len(g),1),C.Counter(r["blind_pci"] for r in good).most_common(8)))
