"""Validate numerical WKB against the exact transmission of the symmetric Eckart
(Poschl-Teller) barrier V(x)=V0*sech^2(a x), V0=Vmax, a=1/(2L).
Exact (Landau & Lifshitz): T = sinh^2(pi k/a) / (sinh^2(pi k/a) + cosh^2(0.5*pi*sqrt(8 m V0/(hbar^2 a^2) - 1)))."""
import numpy as np
from wkb_core import *
rows=[]
for Vmax_eV,LA in ((0.9,0.4),(0.9,0.6),(1.5,0.25),(0.7,0.5)):
    V1=0.0; V2=solve_V2_for_barrier_height(Vmax_eV*EV_TO_J,V1); L=LA*1e-10
    a=1/(2*L); V0=Vmax_eV*EV_TO_J
    for m,name in ((M_PROTON,'H'),(M_DEUTERON,'D')):
        for frac in (0.3,0.5,0.7,0.9):
            E=frac*V0; k=np.sqrt(2*m*E)/HBAR
            ex=np.sinh(np.pi*k/a)**2/(np.sinh(np.pi*k/a)**2+np.cosh(0.5*np.pi*np.sqrt(8*m*V0/(HBAR**2*a**2)-1))**2)
            wk=wkb_transmission(E,V1,V2,L,m)
            rows.append((Vmax_eV,LA,name,frac,ex,wk,wk/ex))
print("Vmax L iso E/V0   exact      WKB      WKB/exact")
for r in rows: print("%.1f %.2f %s %.1f %.3e %.3e %.2f"%r)
import numpy as np
rat=np.array([r[-1] for r in rows]); 
print("median ratio",np.median(rat),"range",rat.min(),rat.max())
