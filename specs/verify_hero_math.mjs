// Regression guard for the WSGA hero math (seed 42, N=160).
// Mirrors docs/index.html: IPW balance (std-diff) AND synthetic effect gap (Delta).
// Run with: node specs/verify_hero_math.mjs
function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
function dnorm(x,m,s){return Math.exp(-0.5*((x-m)/s)**2)/(s*Math.sqrt(2*Math.PI));}
const MU0=0.42,MU1=0.58,SD=0.17,N=160,SEED=42;
const ALPHA=0.1895,BETA=0.0924,GAMMA=-0.3340;
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
const stdDiff=t=>Math.abs(wmean(G1,t)-wmean(G0,t))/Math.sqrt((wsd(G0,t)**2+wsd(G1,t)**2)/2);
const tau=(p,g)=>ALPHA+BETA*g+GAMMA*p.m;
const tauMean=(a,g,t)=>{let s=0,sw=0;for(const p of a){let x=w(p,t);s+=x*tau(p,g);sw+=x;}return s/sw;};
const effectGap=t=>tauMean(G1,1,t)-tauMean(G0,0,t);
const sd0=stdDiff(0),sd1=stdDiff(1),g0=effectGap(0),g1=effectGap(1);
console.log('std-diff  raw',sd0.toFixed(3),'wtd',sd1.toFixed(3),' | effect gap  raw',(g0*100).toFixed(1)+'pp','wtd',(g1*100).toFixed(1)+'pp');
if(sd0<0.90) throw new Error('FAIL std-diff raw ~0.94, got '+sd0.toFixed(3));
if(sd1>0.08) throw new Error('FAIL std-diff wtd ~0.04, got '+sd1.toFixed(3));
if(g0<0.035||g0>0.045) throw new Error('FAIL effect gap raw ~0.040, got '+g0.toFixed(3));
if(g1<0.085||g1>0.095) throw new Error('FAIL effect gap wtd ~0.090, got '+g1.toFixed(3));
console.log('PASS');
