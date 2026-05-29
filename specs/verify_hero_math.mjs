// Regression guard for the WSGA hero IPW math.
// Mirrors the algorithm contract in docs/index.html. Run with: node specs/verify_hero_math.mjs
function rnorm(m,s){let u=0,v=0;while(!u)u=Math.random();while(!v)v=Math.random();return m+s*Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}
function dnorm(x,m,s){return Math.exp(-0.5*((x-m)/s)**2)/(s*Math.sqrt(2*Math.PI));}
const MU0=0.42,MU1=0.58,SD=0.17,N=1000;  // larger N for deterministic threshold; page uses N=160
const G0=[],G1=[];
for(let i=0;i<N;i++)G0.push({m:rnorm(MU0,SD)});
for(let i=0;i<N;i++)G1.push({m:rnorm(MU1,SD)});
for(const [arr,g] of [[G0,0],[G1,1]]){
  for(const p of arr){const f0=dnorm(p.m,MU0,SD),f1=dnorm(p.m,MU1,SD),fp=0.5*(f0+f1);p.tw=fp/(g===0?f0:f1);}
  const mean=arr.reduce((a,p)=>a+p.tw,0)/arr.length; arr.forEach(p=>p.tw/=mean);
}
const w=(p,t)=>1+t*(p.tw-1);
const wmean=(a,t)=>{let s=0,sw=0;for(const p of a){let x=w(p,t);s+=x*p.m;sw+=x;}return s/sw;};
const wsd=(a,t)=>{let mu=wmean(a,t),s=0,sw=0;for(const p of a){let x=w(p,t);s+=x*(p.m-mu)**2;sw+=x;}return Math.sqrt(s/sw);};
const sd=t=>Math.abs(wmean(G1,t)-wmean(G0,t))/Math.sqrt((wsd(G0,t)**2+wsd(G1,t)**2)/2);
const raw=sd(0), wtd=sd(1);
console.log('unweighted std-diff:',raw.toFixed(3),' weighted std-diff:',wtd.toFixed(3));
if(raw < 0.7) throw new Error('FAIL: unweighted std-diff should be large (~0.97), got '+raw.toFixed(3));
if(wtd > 0.20) throw new Error('FAIL: weighted std-diff should collapse (<0.20, ~0.05), got '+wtd.toFixed(3));
console.log('PASS');
