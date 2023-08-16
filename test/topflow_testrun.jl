using SimpleTopOpt
using SimpleTopOpt.TopFlow

Lx = 1.0; Ly = 1.0; nely = 30
volfrac = 1/3 ; Uin = 1e0; rho = 1e0
mu = 1e0; conit = 50

tfdc = TopflowDomain(Lx, Ly, nely)
fea = SimpleTopOpt.TopflowFEA(tfdc)
optimizer = OCParameters(200, 0.2)
dpbc = SimpleTopOpt.DoublePipeBC(tfdc, fea, Uin)
dpc = DoublePipeContainer(tfdc, volfrac, optimizer, Uin, rho, mu)

t_xPhys = SimpleTopOpt.TopFlow.topflow(dpc)