// Regression guard for the WSGA hero IPW math.
// Mirrors the SEEDED algorithm in docs/index.html (seed 42, N=160) so it tests
// the exact instance the page renders. Run with: node specs/verify_hero_math.mjs
function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
function dnorm(x,m,s){return Math.exp(-0.5*((x-m)/s)**2)/(s*Math.sqrt(2*Math.PI));}
const MU0=0.42,MU1=0.58,SD=0.17,N=160,SEED=42;
const rng=mulberry32(SEED);
function rnorm(m,s){let u=0,v=0;while(!u)u=rng();while(!v)v=rng();return m+s*Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}
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
console.log('seed',SEED,'unweighted std-diff:',raw.toFixed(3),' weighted std-diff:',wtd.toFixed(3));
if(raw < 0.90) throw new Error('FAIL: unweighted std-diff should be ~0.94, got '+raw.toFixed(3));
if(wtd > 0.08) throw new Error('FAIL: weighted std-diff should collapse to ~0.04, got '+wtd.toFixed(3));
console.log('PASS');
