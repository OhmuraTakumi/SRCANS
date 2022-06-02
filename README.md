# SRCANS

Special Relativistic Magnetohydridynamics Code : SR-CANS+

----
2022/5/30

Output valiables: ro(rest mass density of fluid),pr,vx,vy,vz(4-velocity.),bx,by,bz(laboratory B-vector)

HLL Rieman Solver, HLLC Rieman Solver(Hokkila and Janhunen 2007, Kim and Balsara 2014 (注：エラーあるかも))

Primitve Ricover (Mignone, Bodo 2006, first guess: Mignone, Mckinney 2007)

Wave speed estimate (Approx. Equation by Leismann et al. 2005) -> 4次方程式を解くように変更すべき

MUSCL(MC, Van leer, minmod), MP5 <- 特性変数補間無し

Shock flattening (Mignone and Bodo 2005): Strong shock zones --> Use minmode limiter

Div. B cleaning, MPI, Cartesian/cylindrical grid,他CANS+を流用
