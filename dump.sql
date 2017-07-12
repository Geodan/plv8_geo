\echo Use "CREATE EXTENSION plv8geo" to load this file. \quit
CREATE TABLE IF NOT EXISTS public.plv8_modules(modname text primary key, load_on_start boolean, code text);
CREATE SCHEMA IF NOT EXISTS plv8;

INSERT INTO plv8_modules (modname, load_on_start, code) VALUES ('d3_force', true, '// https://d3js.org/d3-force/ Version 1.0.6. Copyright 2017 Mike Bostock.
!function(n,t){"object"==typeof exports&&"undefined"!=typeof module?t(exports,require("d3-quadtree"),require("d3-collection"),require("d3-dispatch"),require("d3-timer")):"function"==typeof define&&define.amd?define(["exports","d3-quadtree","d3-collection","d3-dispatch","d3-timer"],t):t(n.d3=n.d3||{},n.d3,n.d3,n.d3,n.d3)}(this,function(n,t,e,r,i){"use strict";function u(n){return n.x+n.vx}function o(n){return n.y+n.vy}function f(n){return n.index}function a(n,t){var e=n.get(t);if(!e)throw new Error("missing: "+t);return e}function c(n){return n.x}function l(n){return n.y}var h=function(n,t){function e(){var e,i,u=r.length,o=0,f=0;for(e=0;e<u;++e)i=r[e],o+=i.x,f+=i.y;for(o=o/u-n,f=f/u-t,e=0;e<u;++e)i=r[e],i.x-=o,i.y-=f}var r;return null==n&&(n=0),null==t&&(t=0),e.initialize=function(n){r=n},e.x=function(t){return arguments.length?(n=+t,e):n},e.y=function(n){return arguments.length?(t=+n,e):t},e},v=function(n){return function(){return n}},d=function(){return 1e-6*(Math.random()-.5)},y=function(n){function e(){function n(n,t,e,r,i){var u=n.data,o=n.r,f=x+o;{if(!u)return t>v+f||r<v-f||e>y+f||i<y-f;if(u.index>h.index){var a=v-u.x-u.vx,l=y-u.y-u.vy,s=a*a+l*l;s<f*f&&(0===a&&(a=d(),s+=a*a),0===l&&(l=d(),s+=l*l),s=(f-(s=Math.sqrt(s)))/s*c,h.vx+=(a*=s)*(f=(o*=o)/(g+o)),h.vy+=(l*=s)*f,u.vx-=a*(f=1-f),u.vy-=l*f)}}}for(var e,i,h,v,y,x,g,s=f.length,p=0;p<l;++p)for(i=t.quadtree(f,u,o).visitAfter(r),e=0;e<s;++e)h=f[e],x=a[h.index],g=x*x,v=h.x+h.vx,y=h.y+h.vy,i.visit(n)}function r(n){if(n.data)return n.r=a[n.data.index];for(var t=n.r=0;t<4;++t)n[t]&&n[t].r>n.r&&(n.r=n[t].r)}function i(){if(f){var t,e,r=f.length;for(a=new Array(r),t=0;t<r;++t)e=f[t],a[e.index]=+n(e,t,f)}}var f,a,c=1,l=1;return"function"!=typeof n&&(n=v(null==n?1:+n)),e.initialize=function(n){f=n,i()},e.iterations=function(n){return arguments.length?(l=+n,e):l},e.strength=function(n){return arguments.length?(c=+n,e):c},e.radius=function(t){return arguments.length?(n="function"==typeof t?t:v(+t),i(),e):n},e},x=function(n){function t(n){return 1/Math.min(y[n.source.index],y[n.target.index])}function r(t){for(var e=0,r=n.length;e<M;++e)for(var i,u,o,f,a,h,v,y=0;y<r;++y)i=n[y],u=i.source,o=i.target,f=o.x+o.vx-u.x-u.vx||d(),a=o.y+o.vy-u.y-u.vy||d(),h=Math.sqrt(f*f+a*a),h=(h-l[y])/h*t*c[y],f*=h,a*=h,o.vx-=f*(v=x[y]),o.vy-=a*v,u.vx+=f*(v=1-v),u.vy+=a*v}function i(){if(h){var t,r,i=h.length,f=n.length,v=e.map(h,g);for(t=0,y=new Array(i);t<f;++t)r=n[t],r.index=t,"object"!=typeof r.source&&(r.source=a(v,r.source)),"object"!=typeof r.target&&(r.target=a(v,r.target)),y[r.source.index]=(y[r.source.index]||0)+1,y[r.target.index]=(y[r.target.index]||0)+1;for(t=0,x=new Array(f);t<f;++t)r=n[t],x[t]=y[r.source.index]/(y[r.source.index]+y[r.target.index]);c=new Array(f),u(),l=new Array(f),o()}}function u(){if(h)for(var t=0,e=n.length;t<e;++t)c[t]=+s(n[t],t,n)}function o(){if(h)for(var t=0,e=n.length;t<e;++t)l[t]=+p(n[t],t,n)}var c,l,h,y,x,g=f,s=t,p=v(30),M=1;return null==n&&(n=[]),r.initialize=function(n){h=n,i()},r.links=function(t){return arguments.length?(n=t,i(),r):n},r.id=function(n){return arguments.length?(g=n,r):g},r.iterations=function(n){return arguments.length?(M=+n,r):M},r.strength=function(n){return arguments.length?(s="function"==typeof n?n:v(+n),u(),r):s},r.distance=function(n){return arguments.length?(p="function"==typeof n?n:v(+n),o(),r):p},r},g=10,s=Math.PI*(3-Math.sqrt(5)),p=function(n){function t(){u(),p.call("tick",a),c<l&&(x.stop(),p.call("end",a))}function u(){var t,e,r=n.length;for(c+=(v-c)*h,y.each(function(n){n(c)}),t=0;t<r;++t)e=n[t],null==e.fx?e.x+=e.vx*=d:(e.x=e.fx,e.vx=0),null==e.fy?e.y+=e.vy*=d:(e.y=e.fy,e.vy=0)}function o(){for(var t,e=0,r=n.length;e<r;++e){if(t=n[e],t.index=e,isNaN(t.x)||isNaN(t.y)){var i=g*Math.sqrt(e),u=e*s;t.x=i*Math.cos(u),t.y=i*Math.sin(u)}(isNaN(t.vx)||isNaN(t.vy))&&(t.vx=t.vy=0)}}function f(t){return t.initialize&&t.initialize(n),t}var a,c=1,l=.001,h=1-Math.pow(l,1/300),v=0,d=.6,y=e.map(),x=i.timer(t),p=r.dispatch("tick","end");return null==n&&(n=[]),o(),a={tick:u,restart:function(){return x.restart(t),a},stop:function(){return x.stop(),a},nodes:function(t){return arguments.length?(n=t,o(),y.each(f),a):n},alpha:function(n){return arguments.length?(c=+n,a):c},alphaMin:function(n){return arguments.length?(l=+n,a):l},alphaDecay:function(n){return arguments.length?(h=+n,a):+h},alphaTarget:function(n){return arguments.length?(v=+n,a):v},velocityDecay:function(n){return arguments.length?(d=1-n,a):1-d},force:function(n,t){return arguments.length>1?(null==t?y.remove(n):y.set(n,f(t)),a):y.get(n)},find:function(t,e,r){var i,u,o,f,a,c=0,l=n.length;for(null==r?r=1/0:r*=r,c=0;c<l;++c)f=n[c],i=t-f.x,u=e-f.y,(o=i*i+u*u)<r&&(a=f,r=o);return a},on:function(n,t){return arguments.length>1?(p.on(n,t),a):p.on(n)}}},M=function(){function n(n){var e,a=u.length,h=t.quadtree(u,c,l).visitAfter(r);for(f=n,e=0;e<a;++e)o=u[e],h.visit(i)}function e(){if(u){var n,t,e=u.length;for(a=new Array(e),n=0;n<e;++n)t=u[n],a[t.index]=+h(t,n,u)}}function r(n){var t,e,r,i,u,o=0;if(n.length){for(r=i=u=0;u<4;++u)(t=n[u])&&(e=t.value)&&(o+=e,r+=e*t.x,i+=e*t.y);n.x=r/o,n.y=i/o}else{t=n,t.x=t.data.x,t.y=t.data.y;do{o+=a[t.data.index]}while(t=t.next)}n.value=o}function i(n,t,e,r){if(!n.value)return!0;var i=n.x-o.x,u=n.y-o.y,c=r-t,l=i*i+u*u;if(c*c/g<l)return l<x&&(0===i&&(i=d(),l+=i*i),0===u&&(u=d(),l+=u*u),l<y&&(l=Math.sqrt(y*l)),o.vx+=i*n.value*f/l,o.vy+=u*n.value*f/l),!0;if(!(n.length||l>=x)){(n.data!==o||n.next)&&(0===i&&(i=d(),l+=i*i),0===u&&(u=d(),l+=u*u),l<y&&(l=Math.sqrt(y*l)));do{n.data!==o&&(c=a[n.data.index]*f/l,o.vx+=i*c,o.vy+=u*c)}while(n=n.next)}}var u,o,f,a,h=v(-30),y=1,x=1/0,g=.81;return n.initialize=function(n){u=n,e()},n.strength=function(t){return arguments.length?(h="function"==typeof t?t:v(+t),e(),n):h},n.distanceMin=function(t){return arguments.length?(y=t*t,n):Math.sqrt(y)},n.distanceMax=function(t){return arguments.length?(x=t*t,n):Math.sqrt(x)},n.theta=function(t){return arguments.length?(g=t*t,n):Math.sqrt(g)},n},q=function(n){function t(n){for(var t,e=0,o=r.length;e<o;++e)t=r[e],t.vx+=(u[e]-t.x)*i[e]*n}function e(){if(r){var t,e=r.length;for(i=new Array(e),u=new Array(e),t=0;t<e;++t)i[t]=isNaN(u[t]=+n(r[t],t,r))?0:+o(r[t],t,r)}}var r,i,u,o=v(.1);return"function"!=typeof n&&(n=v(null==n?0:+n)),t.initialize=function(n){r=n,e()},t.strength=function(n){return arguments.length?(o="function"==typeof n?n:v(+n),e(),t):o},t.x=function(r){return arguments.length?(n="function"==typeof r?r:v(+r),e(),t):n},t},w=function(n){function t(n){for(var t,e=0,o=r.length;e<o;++e)t=r[e],t.vy+=(u[e]-t.y)*i[e]*n}function e(){if(r){var t,e=r.length;for(i=new Array(e),u=new Array(e),t=0;t<e;++t)i[t]=isNaN(u[t]=+n(r[t],t,r))?0:+o(r[t],t,r)}}var r,i,u,o=v(.1);return"function"!=typeof n&&(n=v(null==n?0:+n)),t.initialize=function(n){r=n,e()},t.strength=function(n){return arguments.length?(o="function"==typeof n?n:v(+n),e(),t):o},t.y=function(r){return arguments.length?(n="function"==typeof r?r:v(+r),e(),t):n},t};n.forceCenter=h,n.forceCollide=y,n.forceLink=x,n.forceManyBody=M,n.forceSimulation=p,n.forceX=q,n.forceY=w,Object.defineProperty(n,"__esModule",{value:!0})});');
INSERT INTO plv8_modules (modname, load_on_start, code) VALUES ('d3_geo', true, '// https://d3js.org/d3-geo/ Version 1.6.4. Copyright 2017 Mike Bostock.
(function(n,t){"object"==typeof exports&&"undefined"!=typeof module?t(exports,require("d3-array")):"function"==typeof define&&define.amd?define(["exports","d3-array"],t):t(n.d3=n.d3||{},n.d3)})(this,function(n,t){"use strict";function r(){this.reset()}function i(n,t,r){var i=n.s=t+r,e=i-t,o=i-e;n.t=t-o+(r-e)}function e(n){return n>1?0:n<-1?Dt:Math.acos(n)}function o(n){return n>1?Ut:n<-1?-Ut:Math.asin(n)}function u(n){return(n=rr(n/2))*n}function c(){}function a(n,t){n&&cr.hasOwnProperty(n.type)&&cr[n.type](n,t)}function l(n,t,r){var i,e=-1,o=n.length-r;for(t.lineStart();++e<o;)i=n[e],t.point(i[0],i[1],i[2]);t.lineEnd()}function f(n,t){var r=-1,i=n.length;for(t.polygonStart();++r<i;)l(n[r],t,1);t.polygonEnd()}function s(){sr.point=h}function p(){g(lt,ft)}function h(n,t){sr.point=g,lt=n,ft=t,n*=Ht,t*=Ht,st=n,pt=Vt(t=t/2+Xt),ht=rr(t)}function g(n,t){n*=Ht,t*=Ht,t=t/2+Xt;var r=n-st,i=r>=0?1:-1,e=i*r,o=Vt(t),u=rr(t),c=ht*u,a=pt*o+c*Vt(e),l=c*i*rr(e);lr.add(Qt(l,a)),st=n,pt=o,ht=u}function v(n){return[Qt(n[1],n[0]),o(n[2])]}function d(n){var t=n[0],r=n[1],i=Vt(r);return[i*Vt(t),i*rr(t),rr(r)]}function E(n,t){return n[0]*t[0]+n[1]*t[1]+n[2]*t[2]}function y(n,t){return[n[1]*t[2]-n[2]*t[1],n[2]*t[0]-n[0]*t[2],n[0]*t[1]-n[1]*t[0]]}function S(n,t){n[0]+=t[0],n[1]+=t[1],n[2]+=t[2]}function m(n,t){return[n[0]*t,n[1]*t,n[2]*t]}function M(n){var t=er(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);n[0]/=t,n[1]/=t,n[2]/=t}function x(n,t){xt.push(_t=[gt=n,dt=n]),t<vt&&(vt=t),t>Et&&(Et=t)}function _(n,t){var r=d([n*Ht,t*Ht]);if(Mt){var i=y(Mt,r),e=[i[1],-i[0],0],o=y(e,i);M(o),o=v(o);var u,c=n-yt,a=c>0?1:-1,l=o[0]*Zt*a,f=Jt(c)>180;f^(a*yt<l&&l<a*n)?(u=o[1]*Zt)>Et&&(Et=u):(l=(l+360)%360-180,f^(a*yt<l&&l<a*n)?(u=-o[1]*Zt)<vt&&(vt=u):(t<vt&&(vt=t),t>Et&&(Et=t))),f?n<yt?C(gt,n)>C(gt,dt)&&(dt=n):C(n,dt)>C(gt,dt)&&(gt=n):dt>=gt?(n<gt&&(gt=n),n>dt&&(dt=n)):n>yt?C(gt,n)>C(gt,dt)&&(dt=n):C(n,dt)>C(gt,dt)&&(gt=n)}else xt.push(_t=[gt=n,dt=n]);t<vt&&(vt=t),t>Et&&(Et=t),Mt=r,yt=n}function N(){gr.point=_}function w(){_t[0]=gt,_t[1]=dt,gr.point=x,Mt=null}function P(n,t){if(Mt){var r=n-yt;hr.add(Jt(r)>180?r+(r>0?360:-360):r)}else St=n,mt=t;sr.point(n,t),_(n,t)}function R(){sr.lineStart()}function A(){P(St,mt),sr.lineEnd(),Jt(hr)>Bt&&(gt=-(dt=180)),_t[0]=gt,_t[1]=dt,Mt=null}function C(n,t){return(t-=n)<0?t+360:t}function j(n,t){return n[0]-t[0]}function q(n,t){return n[0]<=n[1]?n[0]<=t&&t<=n[1]:t<n[0]||n[1]<t}function z(n,t){n*=Ht,t*=Ht;var r=Vt(t);b(r*Vt(n),r*rr(n),rr(t))}function b(n,t,r){++Nt,Pt+=(n-Pt)/Nt,Rt+=(t-Rt)/Nt,At+=(r-At)/Nt}function L(){dr.point=O}function O(n,t){n*=Ht,t*=Ht;var r=Vt(t);Tt=r*Vt(n),kt=r*rr(n),Ft=rr(t),dr.point=G,b(Tt,kt,Ft)}function G(n,t){n*=Ht,t*=Ht;var r=Vt(t),i=r*Vt(n),e=r*rr(n),o=rr(t),u=Qt(er((u=kt*o-Ft*e)*u+(u=Ft*i-Tt*o)*u+(u=Tt*e-kt*i)*u),Tt*i+kt*e+Ft*o);wt+=u,Ct+=u*(Tt+(Tt=i)),jt+=u*(kt+(kt=e)),qt+=u*(Ft+(Ft=o)),b(Tt,kt,Ft)}function T(){dr.point=z}function k(){dr.point=I}function F(){B(Ot,Gt),dr.point=z}function I(n,t){Ot=n,Gt=t,n*=Ht,t*=Ht,dr.point=B;var r=Vt(t);Tt=r*Vt(n),kt=r*rr(n),Ft=rr(t),b(Tt,kt,Ft)}function B(n,t){n*=Ht,t*=Ht;var r=Vt(t),i=r*Vt(n),e=r*rr(n),u=rr(t),c=kt*u-Ft*e,a=Ft*i-Tt*u,l=Tt*e-kt*i,f=er(c*c+a*a+l*l),s=o(f),p=f&&-s/f;zt+=p*c,bt+=p*a,Lt+=p*l,wt+=s,Ct+=s*(Tt+(Tt=i)),jt+=s*(kt+(kt=e)),qt+=s*(Ft+(Ft=u)),b(Tt,kt,Ft)}function D(n,t){return[n>Dt?n-Yt:n<-Dt?n+Yt:n,t]}function U(n,t,r){return(n%=Yt)?t||r?Sr(Y(n),Z(t,r)):Y(n):t||r?Z(t,r):D}function X(n){return function(t,r){return t+=n,[t>Dt?t-Yt:t<-Dt?t+Yt:t,r]}}function Y(n){var t=X(n);return t.invert=X(-n),t}function Z(n,t){function r(n,t){var r=Vt(t),a=Vt(n)*r,l=rr(n)*r,f=rr(t),s=f*i+a*e;return[Qt(l*u-s*c,a*i-f*e),o(s*u+l*c)]}var i=Vt(n),e=rr(n),u=Vt(t),c=rr(t);return r.invert=function(n,t){var r=Vt(t),a=Vt(n)*r,l=rr(n)*r,f=rr(t),s=f*u-l*c;return[Qt(l*u+f*c,a*i+s*e),o(s*i-a*e)]},r}function H(n,t,r,i,e,o){if(r){var u=Vt(t),c=rr(t),a=i*r;null==e?(e=t+i*Yt,o=t-a/2):(e=J(u,e),o=J(u,o),(i>0?e<o:e>o)&&(e+=i*Yt));for(var l,f=e;i>0?f>o:f<o;f-=a)l=v([u,-c*Vt(f),-c*rr(f)]),n.point(l[0],l[1])}}function J(n,t){t=d(t),t[0]-=n,M(t);var r=e(-t[1]);return((-t[2]<0?-r:r)+Yt-Bt)%Yt}function K(n,t,r,i){this.x=n,this.z=t,this.o=r,this.e=i,this.v=!1,this.n=this.p=null}function Q(n){if(t=n.length){for(var t,r,i=0,e=n[0];++i<t;)e.n=r=n[i],r.p=e,e=r;e.n=r=n[0],r.p=e}}function V(n,r,i,e){function o(t,o){return n<=t&&t<=i&&r<=o&&o<=e}function u(t,o,u,a){var f=0,s=0;if(null==t||(f=c(t,u))!==(s=c(o,u))||l(t,o)<0^u>0)do{a.point(0===f||3===f?n:i,f>1?e:r)}while((f=(f+u+4)%4)!==s);else a.point(o[0],o[1])}function c(t,e){return Jt(t[0]-n)<Bt?e>0?0:3:Jt(t[0]-i)<Bt?e>0?2:1:Jt(t[1]-r)<Bt?e>0?1:0:e>0?3:2}function a(n,t){return l(n.x,t.x)}function l(n,t){var r=c(n,1),i=c(t,1);return r!==i?r-i:0===r?t[1]-n[1]:1===r?n[0]-t[0]:2===r?n[1]-t[1]:t[0]-n[0]}return function(c){function l(n,t){o(n,t)&&R.point(n,t)}function f(){for(var t=0,r=0,i=E.length;r<i;++r)for(var o,u,c=E[r],a=1,l=c.length,f=c[0],s=f[0],p=f[1];a<l;++a)o=s,u=p,f=c[a],s=f[0],p=f[1],u<=e?p>e&&(s-o)*(e-u)>(p-u)*(n-o)&&++t:p<=e&&(s-o)*(e-u)<(p-u)*(n-o)&&--t;return t}function s(){R=A,d=[],E=[],P=!0}function p(){var n=f(),r=P&&n,i=(d=t.merge(d)).length;(r||i)&&(c.polygonStart(),r&&(c.lineStart(),u(null,null,1,c),c.lineEnd()),i&&Gr(d,a,n,u,c),c.polygonEnd()),R=c,d=E=y=null}function h(){C.point=v,E&&E.push(y=[]),w=!0,N=!1,x=_=NaN}function g(){d&&(v(S,m),M&&N&&A.rejoin(),d.push(A.result())),C.point=l,N&&R.lineEnd()}function v(t,u){var c=o(t,u);if(E&&y.push([t,u]),w)S=t,m=u,M=c,w=!1,c&&(R.lineStart(),R.point(t,u));else if(c&&N)R.point(t,u);else{var a=[x=Math.max(kr,Math.min(Tr,x)),_=Math.max(kr,Math.min(Tr,_))],l=[t=Math.max(kr,Math.min(Tr,t)),u=Math.max(kr,Math.min(Tr,u))];Lr(a,l,n,r,i,e)?(N||(R.lineStart(),R.point(a[0],a[1])),R.point(l[0],l[1]),c||R.lineEnd(),P=!1):c&&(R.lineStart(),R.point(t,u),P=!1)}x=t,_=u,N=c}var d,E,y,S,m,M,x,_,N,w,P,R=c,A=br(),C={point:l,lineStart:h,lineEnd:g,polygonStart:s,polygonEnd:p};return C}}function W(){Ur.point=nn,Ur.lineEnd=$}function $(){Ur.point=Ur.lineEnd=c}function nn(n,t){n*=Ht,t*=Ht,mr=n,Mr=rr(t),xr=Vt(t),Ur.point=tn}function tn(n,t){n*=Ht,t*=Ht;var r=rr(t),i=Vt(t),e=Jt(n-mr),o=Vt(e),u=rr(e),c=i*u,a=xr*r-Mr*i*o,l=Mr*r+xr*i*o;Dr.add(Qt(er(c*c+a*a),l)),mr=n,Mr=r,xr=i}function rn(n,t){return!(!n||!Kr.hasOwnProperty(n.type))&&Kr[n.type](n,t)}function en(n,t){return 0===Hr(n,t)}function on(n,t){var r=Hr(n[0],n[1]);return Hr(n[0],t)+Hr(t,n[1])<=r+Bt}function un(n,t){return!!Br(n.map(cn),an(t))}function cn(n){return n=n.map(an),n.pop(),n}function an(n){return[n[0]*Ht,n[1]*Ht]}function ln(n,r,i){var e=t.range(n,r-Bt,i).concat(r);return function(n){return e.map(function(t){return[n,t]})}}function fn(n,r,i){var e=t.range(n,r-Bt,i).concat(r);return function(n){return e.map(function(t){return[t,n]})}}function sn(){function n(){return{type:"MultiLineString",coordinates:r()}}function r(){return t.range(Wt(u/E)*E,o,E).map(h).concat(t.range(Wt(f/y)*y,l,y).map(g)).concat(t.range(Wt(e/v)*v,i,v).filter(function(n){return Jt(n%E)>Bt}).map(s)).concat(t.range(Wt(a/d)*d,c,d).filter(function(n){return Jt(n%y)>Bt}).map(p))}var i,e,o,u,c,a,l,f,s,p,h,g,v=10,d=v,E=90,y=360,S=2.5;return n.lines=function(){return r().map(function(n){return{type:"LineString",coordinates:n}})},n.outline=function(){return{type:"Polygon",coordinates:[h(u).concat(g(l).slice(1),h(o).reverse().slice(1),g(f).reverse().slice(1))]}},n.extent=function(t){return arguments.length?n.extentMajor(t).extentMinor(t):n.extentMinor()},n.extentMajor=function(t){return arguments.length?(u=+t[0][0],o=+t[1][0],f=+t[0][1],l=+t[1][1],u>o&&(t=u,u=o,o=t),f>l&&(t=f,f=l,l=t),n.precision(S)):[[u,f],[o,l]]},n.extentMinor=function(t){return arguments.length?(e=+t[0][0],i=+t[1][0],a=+t[0][1],c=+t[1][1],e>i&&(t=e,e=i,i=t),a>c&&(t=a,a=c,c=t),n.precision(S)):[[e,a],[i,c]]},n.step=function(t){return arguments.length?n.stepMajor(t).stepMinor(t):n.stepMinor()},n.stepMajor=function(t){return arguments.length?(E=+t[0],y=+t[1],n):[E,y]},n.stepMinor=function(t){return arguments.length?(v=+t[0],d=+t[1],n):[v,d]},n.precision=function(t){return arguments.length?(S=+t,s=ln(a,c,90),p=fn(e,i,S),h=ln(f,l,90),g=fn(u,o,S),n):S},n.extentMajor([[-180,-90+Bt],[180,90-Bt]]).extentMinor([[-180,-80-Bt],[180,80+Bt]])}function pn(){return sn()()}function hn(){ti.point=gn}function gn(n,t){ti.point=vn,_r=wr=n,Nr=Pr=t}function vn(n,t){ni.add(Pr*n-wr*t),wr=n,Pr=t}function dn(){vn(_r,Nr)}function En(n,t){n<ri&&(ri=n),n>ei&&(ei=n),t<ii&&(ii=t),t>oi&&(oi=t)}function yn(n,t){ci+=n,ai+=t,++li}function Sn(){di.point=mn}function mn(n,t){di.point=Mn,yn(Cr=n,jr=t)}function Mn(n,t){var r=n-Cr,i=t-jr,e=er(r*r+i*i);fi+=e*(Cr+n)/2,si+=e*(jr+t)/2,pi+=e,yn(Cr=n,jr=t)}function xn(){di.point=yn}function _n(){di.point=wn}function Nn(){Pn(Rr,Ar)}function wn(n,t){di.point=Pn,yn(Rr=Cr=n,Ar=jr=t)}function Pn(n,t){var r=n-Cr,i=t-jr,e=er(r*r+i*i);fi+=e*(Cr+n)/2,si+=e*(jr+t)/2,pi+=e,e=jr*n-Cr*t,hi+=e*(Cr+n),gi+=e*(jr+t),vi+=3*e,yn(Cr=n,jr=t)}function Rn(n){this._context=n}function An(n,t){_i.point=Cn,yi=mi=n,Si=Mi=t}function Cn(n,t){mi-=n,Mi-=t,xi.add(er(mi*mi+Mi*Mi)),mi=n,Mi=t}function jn(){this._string=[]}function qn(n){return"m0,"+n+"a"+n+","+n+" 0 1,1 0,"+-2*n+"a"+n+","+n+" 0 1,1 0,"+2*n+"z"}function zn(n){return n.length>1}function bn(n,t){return((n=n.x)[0]<0?n[1]-Ut-Bt:Ut-n[1])-((t=t.x)[0]<0?t[1]-Ut-Bt:Ut-t[1])}function Ln(n){var t,r=NaN,i=NaN,e=NaN;return{lineStart:function(){n.lineStart(),t=1},point:function(o,u){var c=o>0?Dt:-Dt,a=Jt(o-r);Jt(a-Dt)<Bt?(n.point(r,i=(i+u)/2>0?Ut:-Ut),n.point(e,i),n.lineEnd(),n.lineStart(),n.point(c,i),n.point(o,i),t=0):e!==c&&a>=Dt&&(Jt(r-e)<Bt&&(r-=e*Bt),Jt(o-c)<Bt&&(o-=c*Bt),i=On(r,i,o,u),n.point(e,i),n.lineEnd(),n.lineStart(),n.point(c,i),t=0),n.point(r=o,i=u),e=c},lineEnd:function(){n.lineEnd(),r=i=NaN},clean:function(){return 2-t}}}function On(n,t,r,i){var e,o,u=rr(n-r);return Jt(u)>Bt?Kt((rr(t)*(o=Vt(i))*rr(r)-rr(i)*(e=Vt(t))*rr(n))/(e*o*u)):(t+i)/2}function Gn(n,t,r,i){var e;if(null==n)e=r*Ut,i.point(-Dt,e),i.point(0,e),i.point(Dt,e),i.point(Dt,0),i.point(Dt,-e),i.point(0,-e),i.point(-Dt,-e),i.point(-Dt,0),i.point(-Dt,e);else if(Jt(n[0]-t[0])>Bt){var o=n[0]<t[0]?Dt:-Dt;e=r*o/2,i.point(-o,e),i.point(0,e),i.point(o,e)}else i.point(t[0],t[1])}function Tn(n){return function(t){var r=new kn;for(var i in n)r[i]=n[i];return r.stream=t,r}}function kn(){}function Fn(n,t,r){var i=t[1][0]-t[0][0],e=t[1][1]-t[0][1],o=n.clipExtent&&n.clipExtent();n.scale(150).translate([0,0]),null!=o&&n.clipExtent(null),ar(r,n.stream(ui));var u=ui.result(),c=Math.min(i/(u[1][0]-u[0][0]),e/(u[1][1]-u[0][1])),a=+t[0][0]+(i-c*(u[1][0]+u[0][0]))/2,l=+t[0][1]+(e-c*(u[1][1]+u[0][1]))/2;return null!=o&&n.clipExtent(o),n.scale(150*c).translate([a,l])}function In(n,t,r){return Fn(n,[[0,0],t],r)}function Bn(n){return Tn({point:function(t,r){t=n(t,r),this.stream.point(t[0],t[1])}})}function Dn(n,t){function r(i,e,u,c,a,l,f,s,p,h,g,v,d,E){var y=f-i,S=s-e,m=y*y+S*S;if(m>4*t&&d--){var M=c+h,x=a+g,_=l+v,N=er(M*M+x*x+_*_),w=o(_/=N),P=Jt(Jt(_)-1)<Bt||Jt(u-p)<Bt?(u+p)/2:Qt(x,M),R=n(P,w),A=R[0],C=R[1],j=A-i,q=C-e,z=S*j-y*q;(z*z/m>t||Jt((y*j+S*q)/m-.5)>.3||c*h+a*g+l*v<ji)&&(r(i,e,u,c,a,l,A,C,P,M/=N,x/=N,_,d,E),E.point(A,C),r(A,C,P,M,x,_,f,s,p,h,g,v,d,E))}}return function(t){function i(r,i){r=n(r,i),t.point(r[0],r[1])}function e(){y=NaN,_.point=o,t.lineStart()}function o(i,e){var o=d([i,e]),u=n(i,e);r(y,S,E,m,M,x,y=u[0],S=u[1],E=i,m=o[0],M=o[1],x=o[2],Ci,t),t.point(y,S)}function u(){_.point=i,t.lineEnd()}function c(){e(),_.point=a,_.lineEnd=l}function a(n,t){o(f=n,t),s=y,p=S,h=m,g=M,v=x,_.point=o}function l(){r(y,S,E,m,M,x,s,p,f,h,g,v,Ci,t),_.lineEnd=u,u()}var f,s,p,h,g,v,E,y,S,m,M,x,_={point:i,lineStart:e,lineEnd:u,polygonStart:function(){t.polygonStart(),_.lineStart=c},polygonEnd:function(){t.polygonEnd(),_.lineStart=e}};return _}}function Un(n){return Xn(function(){return n})()}function Xn(n){function t(n){return n=f(n[0]*Ht,n[1]*Ht),[n[0]*d+c,a-n[1]*d]}function r(n){return(n=f.invert((n[0]-c)/d,(a-n[1])/d))&&[n[0]*Zt,n[1]*Zt]}function i(n,t){return n=u(n,t),[n[0]*d+c,a-n[1]*d]}function e(){f=Sr(l=U(M,x,_),u);var n=u(S,m);return c=E-n[0]*d,a=y+n[1]*d,o()}function o(){return g=v=null,t}var u,c,a,l,f,s,p,h,g,v,d=150,E=480,y=250,S=0,m=0,M=0,x=0,_=0,N=null,w=Pi,P=null,R=Wr,A=.5,C=qi(i,A);return t.stream=function(n){return g&&v===n?g:g=zi(w(l,C(R(v=n))))},t.clipAngle=function(n){return arguments.length?(w=+n?Ri(N=n*Ht,6*Ht):(N=null,Pi),o()):N*Zt},t.clipExtent=function(n){return arguments.length?(R=null==n?(P=s=p=h=null,Wr):V(P=+n[0][0],s=+n[0][1],p=+n[1][0],h=+n[1][1]),o()):null==P?null:[[P,s],[p,h]]},t.scale=function(n){return arguments.length?(d=+n,e()):d},t.translate=function(n){return arguments.length?(E=+n[0],y=+n[1],e()):[E,y]},t.center=function(n){return arguments.length?(S=n[0]%360*Ht,m=n[1]%360*Ht,e()):[S*Zt,m*Zt]},t.rotate=function(n){return arguments.length?(M=n[0]%360*Ht,x=n[1]%360*Ht,_=n.length>2?n[2]%360*Ht:0,e()):[M*Zt,x*Zt,_*Zt]},t.precision=function(n){return arguments.length?(C=qi(i,A=n*n),o()):er(A)},t.fitExtent=function(n,r){return Fn(t,n,r)},t.fitSize=function(n,r){return In(t,n,r)},function(){return u=n.apply(this,arguments),t.invert=u.invert&&r,e()}}function Yn(n){var t=0,r=Dt/3,i=Xn(n),e=i(t,r);return e.parallels=function(n){return arguments.length?i(t=n[0]*Ht,r=n[1]*Ht):[t*Zt,r*Zt]},e}function Zn(n){function t(n,t){return[n*r,rr(t)/r]}var r=Vt(n);return t.invert=function(n,t){return[n/r,o(t*r)]},t}function Hn(n,t){function r(n,t){var r=er(u-2*e*rr(t))/e;return[r*rr(n*=e),c-r*Vt(n)]}var i=rr(n),e=(i+rr(t))/2;if(Jt(e)<Bt)return Zn(n);var u=1+i*(2*e-i),c=er(u)/e;return r.invert=function(n,t){var r=c-t;return[Qt(n,Jt(r))/e*ir(r),o((u-(n*n+r*r)*e*e)/(2*e))]},r}function Jn(n){var t=n.length;return{point:function(r,i){for(var e=-1;++e<t;)n[e].point(r,i)},sphere:function(){for(var r=-1;++r<t;)n[r].sphere()},lineStart:function(){for(var r=-1;++r<t;)n[r].lineStart()},lineEnd:function(){for(var r=-1;++r<t;)n[r].lineEnd()},polygonStart:function(){for(var r=-1;++r<t;)n[r].polygonStart()},polygonEnd:function(){for(var r=-1;++r<t;)n[r].polygonEnd()}}}function Kn(n){return function(t,r){var i=Vt(t),e=Vt(r),o=n(i*e);return[o*e*rr(t),o*rr(r)]}}function Qn(n){return function(t,r){var i=er(t*t+r*r),e=n(i),u=rr(e),c=Vt(e);return[Qt(t*u,i*c),o(i&&r*u/i)]}}function Vn(n,t){return[n,nr(or((Ut+t)/2))]}function Wn(n){function t(){var t=Dt*c(),u=o(qr(o.rotate()).invert([0,0]));return l(null==f?[[u[0]-t,u[1]-t],[u[0]+t,u[1]+t]]:n===Vn?[[Math.max(u[0]-t,f),r],[Math.min(u[0]+t,i),e]]:[[f,Math.max(u[1]-t,r)],[i,Math.min(u[1]+t,e)]])}var r,i,e,o=Un(n),u=o.center,c=o.scale,a=o.translate,l=o.clipExtent,f=null;return o.scale=function(n){return arguments.length?(c(n),t()):c()},o.translate=function(n){return arguments.length?(a(n),t()):a()},o.center=function(n){return arguments.length?(u(n),t()):u()},o.clipExtent=function(n){return arguments.length?(null==n?f=r=i=e=null:(f=+n[0][0],r=+n[0][1],i=+n[1][0],e=+n[1][1]),t()):null==f?null:[[f,r],[i,e]]},t()}function $n(n){return or((Ut+n)/2)}function nt(n,t){function r(n,t){o>0?t<-Ut+Bt&&(t=-Ut+Bt):t>Ut-Bt&&(t=Ut-Bt);var r=o/tr($n(t),e);return[r*rr(e*n),o-r*Vt(e*n)]}var i=Vt(n),e=n===t?rr(n):nr(i/Vt(t))/nr($n(t)/$n(n)),o=i*tr($n(n),e)/e;return e?(r.invert=function(n,t){var r=o-t,i=ir(e)*er(n*n+r*r);return[Qt(n,Jt(r))/e*ir(r),2*Kt(tr(o/i,1/e))-Ut]},r):Vn}function tt(n,t){return[n,t]}function rt(n,t){function r(n,t){var r=o-t,i=e*n;return[r*rr(i),o-r*Vt(i)]}var i=Vt(n),e=n===t?rr(n):(i-Vt(t))/(t-n),o=i/e+n;return Jt(e)<Bt?tt:(r.invert=function(n,t){var r=o-t;return[Qt(n,Jt(r))/e*ir(r),o-ir(e)*er(n*n+r*r)]},r)}function it(n,t){var r=Vt(t),i=Vt(n)*r;return[r*rr(n)/i,rr(t)/i]}function et(n,t,r,i){return 1===n&&1===t&&0===r&&0===i?Wr:Tn({point:function(e,o){this.stream.point(e*n+r,o*t+i)}})}function ot(n,t){return[Vt(t)*rr(n),rr(t)]}function ut(n,t){var r=Vt(t),i=1+Vt(n)*r;return[r*rr(n)/i,rr(t)/i]}function ct(n,t){return[nr(or((Ut+t)/2)),-n]}var at=function(){return new r};r.prototype={constructor:r,reset:function(){this.s=this.t=0},add:function(n){i(It,n,this.t),i(this,It.s,this.s),this.s?this.t+=It.t:this.s=It.t},valueOf:function(){return this.s}};var lt,ft,st,pt,ht,gt,vt,dt,Et,yt,St,mt,Mt,xt,_t,Nt,wt,Pt,Rt,At,Ct,jt,qt,zt,bt,Lt,Ot,Gt,Tt,kt,Ft,It=new r,Bt=1e-6,Dt=Math.PI,Ut=Dt/2,Xt=Dt/4,Yt=2*Dt,Zt=180/Dt,Ht=Dt/180,Jt=Math.abs,Kt=Math.atan,Qt=Math.atan2,Vt=Math.cos,Wt=Math.ceil,$t=Math.exp,nr=Math.log,tr=Math.pow,rr=Math.sin,ir=Math.sign||function(n){return n>0?1:n<0?-1:0},er=Math.sqrt,or=Math.tan,ur={Feature:function(n,t){a(n.geometry,t)},FeatureCollection:function(n,t){for(var r=n.features,i=-1,e=r.length;++i<e;)a(r[i].geometry,t)}},cr={Sphere:function(n,t){t.sphere()},Point:function(n,t){n=n.coordinates,t.point(n[0],n[1],n[2])},MultiPoint:function(n,t){for(var r=n.coordinates,i=-1,e=r.length;++i<e;)n=r[i],t.point(n[0],n[1],n[2])},LineString:function(n,t){l(n.coordinates,t,0)},MultiLineString:function(n,t){for(var r=n.coordinates,i=-1,e=r.length;++i<e;)l(r[i],t,0)},Polygon:function(n,t){f(n.coordinates,t)},MultiPolygon:function(n,t){for(var r=n.coordinates,i=-1,e=r.length;++i<e;)f(r[i],t)},GeometryCollection:function(n,t){for(var r=n.geometries,i=-1,e=r.length;++i<e;)a(r[i],t)}},ar=function(n,t){n&&ur.hasOwnProperty(n.type)?ur[n.type](n,t):a(n,t)},lr=at(),fr=at(),sr={point:c,lineStart:c,lineEnd:c,polygonStart:function(){lr.reset(),sr.lineStart=s,sr.lineEnd=p},polygonEnd:function(){var n=+lr;fr.add(n<0?Yt+n:n),this.lineStart=this.lineEnd=this.point=c},sphere:function(){fr.add(Yt)}},pr=function(n){return fr.reset(),ar(n,sr),2*fr},hr=at(),gr={point:x,lineStart:N,lineEnd:w,polygonStart:function(){gr.point=P,gr.lineStart=R,gr.lineEnd=A,hr.reset(),sr.polygonStart()},polygonEnd:function(){sr.polygonEnd(),gr.point=x,gr.lineStart=N,gr.lineEnd=w,lr<0?(gt=-(dt=180),vt=-(Et=90)):hr>Bt?Et=90:hr<-Bt&&(vt=-90),_t[0]=gt,_t[1]=dt}},vr=function(n){var t,r,i,e,o,u,c;if(Et=dt=-(gt=vt=1/0),xt=[],ar(n,gr),r=xt.length){for(xt.sort(j),t=1,i=xt[0],o=[i];t<r;++t)e=xt[t],q(i,e[0])||q(i,e[1])?(C(i[0],e[1])>C(i[0],i[1])&&(i[1]=e[1]),C(e[0],i[1])>C(i[0],i[1])&&(i[0]=e[0])):o.push(i=e);for(u=-1/0,r=o.length-1,t=0,i=o[r];t<=r;i=e,++t)e=o[t],(c=C(i[1],e[0]))>u&&(u=c,gt=e[0],dt=i[1])}return xt=_t=null,gt===1/0||vt===1/0?[[NaN,NaN],[NaN,NaN]]:[[gt,vt],[dt,Et]]},dr={sphere:c,point:z,lineStart:L,lineEnd:T,polygonStart:function(){dr.lineStart=k,dr.lineEnd=F},polygonEnd:function(){dr.lineStart=L,dr.lineEnd=T}},Er=function(n){Nt=wt=Pt=Rt=At=Ct=jt=qt=zt=bt=Lt=0,ar(n,dr);var t=zt,r=bt,i=Lt,e=t*t+r*r+i*i;return e<1e-12&&(t=Ct,r=jt,i=qt,wt<Bt&&(t=Pt,r=Rt,i=At),(e=t*t+r*r+i*i)<1e-12)?[NaN,NaN]:[Qt(r,t)*Zt,o(i/er(e))*Zt]},yr=function(n){return function(){return n}},Sr=function(n,t){function r(r,i){return r=n(r,i),t(r[0],r[1])}return n.invert&&t.invert&&(r.invert=function(r,i){return(r=t.invert(r,i))&&n.invert(r[0],r[1])}),r};D.invert=D;var mr,Mr,xr,_r,Nr,wr,Pr,Rr,Ar,Cr,jr,qr=function(n){function t(t){return t=n(t[0]*Ht,t[1]*Ht),t[0]*=Zt,t[1]*=Zt,t}return n=U(n[0]*Ht,n[1]*Ht,n.length>2?n[2]*Ht:0),t.invert=function(t){return t=n.invert(t[0]*Ht,t[1]*Ht),t[0]*=Zt,t[1]*=Zt,t},t},zr=function(){function n(n,t){r.push(n=i(n,t)),n[0]*=Zt,n[1]*=Zt}function t(){var n=e.apply(this,arguments),t=o.apply(this,arguments)*Ht,a=u.apply(this,arguments)*Ht;return r=[],i=U(-n[0]*Ht,-n[1]*Ht,0).invert,H(c,t,a,1),n={type:"Polygon",coordinates:[r]},r=i=null,n}var r,i,e=yr([0,0]),o=yr(90),u=yr(6),c={point:n};return t.center=function(n){return arguments.length?(e="function"==typeof n?n:yr([+n[0],+n[1]]),t):e},t.radius=function(n){return arguments.length?(o="function"==typeof n?n:yr(+n),t):o},t.precision=function(n){return arguments.length?(u="function"==typeof n?n:yr(+n),t):u},t},br=function(){var n,t=[];return{point:function(t,r){n.push([t,r])},lineStart:function(){t.push(n=[])},lineEnd:c,rejoin:function(){t.length>1&&t.push(t.pop().concat(t.shift()))},result:function(){var r=t;return t=[],n=null,r}}},Lr=function(n,t,r,i,e,o){var u,c=n[0],a=n[1],l=t[0],f=t[1],s=0,p=1,h=l-c,g=f-a;if(u=r-c,h||!(u>0)){if(u/=h,h<0){if(u<s)return;u<p&&(p=u)}else if(h>0){if(u>p)return;u>s&&(s=u)}if(u=e-c,h||!(u<0)){if(u/=h,h<0){if(u>p)return;u>s&&(s=u)}else if(h>0){if(u<s)return;u<p&&(p=u)}if(u=i-a,g||!(u>0)){if(u/=g,g<0){if(u<s)return;u<p&&(p=u)}else if(g>0){if(u>p)return;u>s&&(s=u)}if(u=o-a,g||!(u<0)){if(u/=g,g<0){if(u>p)return;u>s&&(s=u)}else if(g>0){if(u<s)return;u<p&&(p=u)}return s>0&&(n[0]=c+s*h,n[1]=a+s*g),p<1&&(t[0]=c+p*h,t[1]=a+p*g),!0}}}}},Or=function(n,t){return Jt(n[0]-t[0])<Bt&&Jt(n[1]-t[1])<Bt},Gr=function(n,t,r,i,e){var o,u,c=[],a=[];if(n.forEach(function(n){if(!((t=n.length-1)<=0)){var t,r,i=n[0],u=n[t];if(Or(i,u)){for(e.lineStart(),o=0;o<t;++o)e.point((i=n[o])[0],i[1]);return void e.lineEnd()}c.push(r=new K(i,n,null,!0)),a.push(r.o=new K(i,null,r,!1)),c.push(r=new K(u,n,null,!1)),a.push(r.o=new K(u,null,r,!0))}}),c.length){for(a.sort(t),Q(c),Q(a),o=0,u=a.length;o<u;++o)a[o].e=r=!r;for(var l,f,s=c[0];;){for(var p=s,h=!0;p.v;)if((p=p.n)===s)return;l=p.z,e.lineStart();do{if(p.v=p.o.v=!0,p.e){if(h)for(o=0,u=l.length;o<u;++o)e.point((f=l[o])[0],f[1]);else i(p.x,p.n.x,1,e);p=p.n}else{if(h)for(l=p.p.z,o=l.length-1;o>=0;--o)e.point((f=l[o])[0],f[1]);else i(p.x,p.p.x,-1,e);p=p.p}p=p.o,l=p.z,h=!h}while(!p.v);e.lineEnd()}}},Tr=1e9,kr=-Tr,Fr=function(){var n,t,r,i=0,e=0,o=960,u=500;return r={stream:function(r){return n&&t===r?n:n=V(i,e,o,u)(t=r)},extent:function(c){return arguments.length?(i=+c[0][0],e=+c[0][1],o=+c[1][0],u=+c[1][1],n=t=null,r):[[i,e],[o,u]]}}},Ir=at(),Br=function(n,t){var r=t[0],i=t[1],e=[rr(r),-Vt(r),0],u=0,c=0;Ir.reset();for(var a=0,l=n.length;a<l;++a)if(s=(f=n[a]).length)for(var f,s,p=f[s-1],h=p[0],g=p[1]/2+Xt,v=rr(g),E=Vt(g),S=0;S<s;++S,h=x,v=N,E=w,p=m){var m=f[S],x=m[0],_=m[1]/2+Xt,N=rr(_),w=Vt(_),P=x-h,R=P>=0?1:-1,A=R*P,C=A>Dt,j=v*N;if(Ir.add(Qt(j*R*rr(A),E*w+j*Vt(A))),u+=C?P+R*Yt:P,C^h>=r^x>=r){var q=y(d(p),d(m));M(q);var z=y(e,q);M(z);var b=(C^P>=0?-1:1)*o(z[2]);(i>b||i===b&&(q[0]||q[1]))&&(c+=C^P>=0?1:-1)}}return(u<-Bt||u<Bt&&Ir<-Bt)^1&c},Dr=at(),Ur={sphere:c,point:c,lineStart:W,lineEnd:c,polygonStart:c,polygonEnd:c},Xr=function(n){return Dr.reset(),ar(n,Ur),+Dr},Yr=[null,null],Zr={type:"LineString",coordinates:Yr},Hr=function(n,t){return Yr[0]=n,Yr[1]=t,Xr(Zr)},Jr={Feature:function(n,t){return rn(n.geometry,t)},FeatureCollection:function(n,t){for(var r=n.features,i=-1,e=r.length;++i<e;)if(rn(r[i].geometry,t))return!0;return!1}},Kr={Sphere:function(){return!0},Point:function(n,t){return en(n.coordinates,t)},MultiPoint:function(n,t){for(var r=n.coordinates,i=-1,e=r.length;++i<e;)if(en(r[i],t))return!0;return!1},LineString:function(n,t){return on(n.coordinates,t)},MultiLineString:function(n,t){for(var r=n.coordinates,i=-1,e=r.length;++i<e;)if(on(r[i],t))return!0;return!1},Polygon:function(n,t){return un(n.coordinates,t)},MultiPolygon:function(n,t){for(var r=n.coordinates,i=-1,e=r.length;++i<e;)if(un(r[i],t))return!0;return!1},GeometryCollection:function(n,t){for(var r=n.geometries,i=-1,e=r.length;++i<e;)if(rn(r[i],t))return!0;return!1}},Qr=function(n,t){return(n&&Jr.hasOwnProperty(n.type)?Jr[n.type]:rn)(n,t)},Vr=function(n,t){var r=n[0]*Ht,i=n[1]*Ht,e=t[0]*Ht,c=t[1]*Ht,a=Vt(i),l=rr(i),f=Vt(c),s=rr(c),p=a*Vt(r),h=a*rr(r),g=f*Vt(e),v=f*rr(e),d=2*o(er(u(c-i)+a*f*u(e-r))),E=rr(d),y=d?function(n){var t=rr(n*=d)/E,r=rr(d-n)/E,i=r*p+t*g,e=r*h+t*v,o=r*l+t*s;return[Qt(e,i)*Zt,Qt(o,er(i*i+e*e))*Zt]}:function(){return[r*Zt,i*Zt]};return y.distance=d,y},Wr=function(n){return n},$r=at(),ni=at(),ti={point:c,lineStart:c,lineEnd:c,polygonStart:function(){ti.lineStart=hn,ti.lineEnd=dn},polygonEnd:function(){ti.lineStart=ti.lineEnd=ti.point=c,$r.add(Jt(ni)),ni.reset()},result:function(){var n=$r/2;return $r.reset(),n}},ri=1/0,ii=ri,ei=-ri,oi=ei,ui={point:En,lineStart:c,lineEnd:c,polygonStart:c,polygonEnd:c,result:function(){var n=[[ri,ii],[ei,oi]];return ei=oi=-(ii=ri=1/0),n}},ci=0,ai=0,li=0,fi=0,si=0,pi=0,hi=0,gi=0,vi=0,di={point:yn,lineStart:Sn,lineEnd:xn,polygonStart:function(){di.lineStart=_n,di.lineEnd=Nn},polygonEnd:function(){di.point=yn,di.lineStart=Sn,di.lineEnd=xn},result:function(){var n=vi?[hi/vi,gi/vi]:pi?[fi/pi,si/pi]:li?[ci/li,ai/li]:[NaN,NaN];return ci=ai=li=fi=si=pi=hi=gi=vi=0,n}};Rn.prototype={_radius:4.5,pointRadius:function(n){return this._radius=n,this},polygonStart:function(){this._line=0},polygonEnd:function(){this._line=NaN},lineStart:function(){this._point=0},lineEnd:function(){0===this._line&&this._context.closePath(),this._point=NaN},point:function(n,t){switch(this._point){case 0:this._context.moveTo(n,t),this._point=1;break;case 1:this._context.lineTo(n,t);break;default:this._context.moveTo(n+this._radius,t),this._context.arc(n,t,this._radius,0,Yt)}},result:c};var Ei,yi,Si,mi,Mi,xi=at(),_i={point:c,lineStart:function(){_i.point=An},lineEnd:function(){Ei&&Cn(yi,Si),_i.point=c},polygonStart:function(){Ei=!0},polygonEnd:function(){Ei=null},result:function(){var n=+xi;return xi.reset(),n}};jn.prototype={_radius:4.5,_circle:qn(4.5),pointRadius:function(n){return(n=+n)!==this._radius&&(this._radius=n,this._circle=null),this},polygonStart:function(){this._line=0},polygonEnd:function(){this._line=NaN},lineStart:function(){this._point=0},lineEnd:function(){0===this._line&&this._string.push("Z"),this._point=NaN},point:function(n,t){switch(this._point){case 0:this._string.push("M",n,",",t),this._point=1;break;case 1:this._string.push("L",n,",",t);break;default:null==this._circle&&(this._circle=qn(this._radius)),this._string.push("M",n,",",t,this._circle)}},result:function(){if(this._string.length){var n=this._string.join("");return this._string=[],n}return null}};var Ni=function(n,t){function r(n){return n&&("function"==typeof o&&e.pointRadius(+o.apply(this,arguments)),ar(n,i(e))),e.result()}var i,e,o=4.5;return r.area=function(n){return ar(n,i(ti)),ti.result()},r.measure=function(n){return ar(n,i(_i)),_i.result()},r.bounds=function(n){return ar(n,i(ui)),ui.result()},r.centroid=function(n){return ar(n,i(di)),di.result()},r.projection=function(t){return arguments.length?(i=null==t?(n=null,Wr):(n=t).stream,r):n},r.context=function(n){return arguments.length?(e=null==n?(t=null,new jn):new Rn(t=n),"function"!=typeof o&&e.pointRadius(o),r):t},r.pointRadius=function(n){return arguments.length?(o="function"==typeof n?n:(e.pointRadius(+n),+n),r):o},r.projection(n).context(t)},wi=function(n,r,i,e){return function(o,u){function c(t,r){var i=o(t,r);n(t=i[0],r=i[1])&&u.point(t,r)}function a(n,t){var r=o(n,t);E.point(r[0],r[1])}function l(){x.point=a,E.lineStart()}function f(){x.point=c,E.lineEnd()}function s(n,t){d.push([n,t]);var r=o(n,t);m.point(r[0],r[1])}function p(){m.lineStart(),d=[]}function h(){s(d[0][0],d[0][1]),m.lineEnd();var n,t,r,i,e=m.clean(),o=S.result(),c=o.length;if(d.pop(),g.push(d),d=null,c)if(1&e){if(r=o[0],(t=r.length-1)>0){for(M||(u.polygonStart(),M=!0),u.lineStart(),n=0;n<t;++n)u.point((i=r[n])[0],i[1]);u.lineEnd()}}else c>1&&2&e&&o.push(o.pop().concat(o.shift())),v.push(o.filter(zn))}var g,v,d,E=r(u),y=o.invert(e[0],e[1]),S=br(),m=r(S),M=!1,x={point:c,lineStart:l,lineEnd:f,polygonStart:function(){x.point=s,x.lineStart=p,x.lineEnd=h,v=[],g=[]},polygonEnd:function(){x.point=c,x.lineStart=l,x.lineEnd=f,v=t.merge(v);var n=Br(g,y);v.length?(M||(u.polygonStart(),M=!0),Gr(v,bn,n,i,u)):n&&(M||(u.polygonStart(),M=!0),u.lineStart(),i(null,null,1,u),u.lineEnd()),M&&(u.polygonEnd(),M=!1),v=g=null},sphere:function(){u.polygonStart(),u.lineStart(),i(null,null,1,u),u.lineEnd(),u.polygonEnd()}};return x}},Pi=wi(function(){return!0},Ln,Gn,[-Dt,-Ut]),Ri=function(n,t){function r(r,i,e,o){H(o,n,t,e,r,i)}function i(n,t){return Vt(n)*Vt(t)>c}function e(n){var t,r,e,c,f;return{lineStart:function(){c=e=!1,f=1},point:function(s,p){var h,g=[s,p],v=i(s,p),d=a?v?0:u(s,p):v?u(s+(s<0?Dt:-Dt),p):0;if(!t&&(c=e=v)&&n.lineStart(),v!==e&&(!(h=o(t,g))||Or(t,h)||Or(g,h))&&(g[0]+=Bt,g[1]+=Bt,v=i(g[0],g[1])),v!==e)f=0,v?(n.lineStart(),h=o(g,t),n.point(h[0],h[1])):(h=o(t,g),n.point(h[0],h[1]),n.lineEnd()),t=h;else if(l&&t&&a^v){var E;d&r||!(E=o(g,t,!0))||(f=0,a?(n.lineStart(),n.point(E[0][0],E[0][1]),n.point(E[1][0],E[1][1]),n.lineEnd()):(n.point(E[1][0],E[1][1]),n.lineEnd(),n.lineStart(),n.point(E[0][0],E[0][1])))}!v||t&&Or(t,g)||n.point(g[0],g[1]),t=g,e=v,r=d},lineEnd:function(){e&&n.lineEnd(),t=null},clean:function(){return f|(c&&e)<<1}}}function o(n,t,r){var i=d(n),e=d(t),o=[1,0,0],u=y(i,e),a=E(u,u),l=u[0],f=a-l*l;if(!f)return!r&&n;var s=c*a/f,p=-c*l/f,h=y(o,u),g=m(o,s);S(g,m(u,p));var M=h,x=E(g,M),_=E(M,M),N=x*x-_*(E(g,g)-1);if(!(N<0)){var w=er(N),P=m(M,(-x-w)/_);if(S(P,g),P=v(P),!r)return P;var R,A=n[0],C=t[0],j=n[1],q=t[1];C<A&&(R=A,A=C,C=R);var z=C-A,b=Jt(z-Dt)<Bt,L=b||z<Bt;if(!b&&q<j&&(R=j,j=q,q=R),L?b?j+q>0^P[1]<(Jt(P[0]-A)<Bt?j:q):j<=P[1]&&P[1]<=q:z>Dt^(A<=P[0]&&P[0]<=C)){var O=m(M,(-x+w)/_);return S(O,g),[P,v(O)]}}}function u(t,r){var i=a?n:Dt-n,e=0;return t<-i?e|=1:t>i&&(e|=2),r<-i?e|=4:r>i&&(e|=8),e}var c=Vt(n),a=c>0,l=Jt(c)>Bt;return wi(i,e,r,a?[0,-n]:[-Dt,n-Dt])},Ai=function(n){return{stream:Tn(n)}};kn.prototype={constructor:kn,point:function(n,t){this.stream.point(n,t)},sphere:function(){this.stream.sphere()},lineStart:function(){this.stream.lineStart()},lineEnd:function(){this.stream.lineEnd()},polygonStart:function(){this.stream.polygonStart()},polygonEnd:function(){this.stream.polygonEnd()}};var Ci=16,ji=Vt(30*Ht),qi=function(n,t){return+t?Dn(n,t):Bn(n)},zi=Tn({point:function(n,t){this.stream.point(n*Ht,t*Ht)}}),bi=function(){return Yn(Hn).scale(155.424).center([0,33.6442])},Li=function(){return bi().parallels([29.5,45.5]).scale(1070).translate([480,250]).rotate([96,0]).center([-.6,38.7])},Oi=function(){function n(n){var t=n[0],r=n[1];return c=null,e.point(t,r),c||(o.point(t,r),c)||(u.point(t,r),c)}function t(){return r=i=null,n}var r,i,e,o,u,c,a=Li(),l=bi().rotate([154,0]).center([-2,58.5]).parallels([55,65]),f=bi().rotate([157,0]).center([-3,19.9]).parallels([8,18]),s={point:function(n,t){c=[n,t]}};return n.invert=function(n){var t=a.scale(),r=a.translate(),i=(n[0]-r[0])/t,e=(n[1]-r[1])/t;return(e>=.12&&e<.234&&i>=-.425&&i<-.214?l:e>=.166&&e<.234&&i>=-.214&&i<-.115?f:a).invert(n)},n.stream=function(n){return r&&i===n?r:r=Jn([a.stream(i=n),l.stream(n),f.stream(n)])},n.precision=function(n){return arguments.length?(a.precision(n),l.precision(n),f.precision(n),t()):a.precision()},n.scale=function(t){return arguments.length?(a.scale(t),l.scale(.35*t),f.scale(t),n.translate(a.translate())):a.scale()},n.translate=function(n){if(!arguments.length)return a.translate();var r=a.scale(),i=+n[0],c=+n[1];return e=a.translate(n).clipExtent([[i-.455*r,c-.238*r],[i+.455*r,c+.238*r]]).stream(s),o=l.translate([i-.307*r,c+.201*r]).clipExtent([[i-.425*r+Bt,c+.12*r+Bt],[i-.214*r-Bt,c+.234*r-Bt]]).stream(s),u=f.translate([i-.205*r,c+.212*r]).clipExtent([[i-.214*r+Bt,c+.166*r+Bt],[i-.115*r-Bt,c+.234*r-Bt]]).stream(s),t()},n.fitExtent=function(t,r){return Fn(n,t,r)},n.fitSize=function(t,r){return In(n,t,r)},n.scale(1070)},Gi=Kn(function(n){return er(2/(1+n))});Gi.invert=Qn(function(n){return 2*o(n/2)});var Ti=function(){return Un(Gi).scale(124.75).clipAngle(179.999)},ki=Kn(function(n){return(n=e(n))&&n/rr(n)});ki.invert=Qn(function(n){return n});var Fi=function(){return Un(ki).scale(79.4188).clipAngle(179.999)};Vn.invert=function(n,t){return[n,2*Kt($t(t))-Ut]};var Ii=function(){return Wn(Vn).scale(961/Yt)},Bi=function(){return Yn(nt).scale(109.5).parallels([30,30])};tt.invert=tt;var Di=function(){return Un(tt).scale(152.63)},Ui=function(){return Yn(rt).scale(131.154).center([0,13.9389])};it.invert=Qn(Kt);var Xi=function(){return Un(it).scale(144.049).clipAngle(60)},Yi=function(){function n(){return e=o=null,u}var t,r,i,e,o,u,c=1,a=0,l=0,f=1,s=1,p=Wr,h=null,g=Wr;return u={stream:function(n){return e&&o===n?e:e=p(g(o=n))},clipExtent:function(e){return arguments.length?(g=null==e?(h=t=r=i=null,Wr):V(h=+e[0][0],t=+e[0][1],r=+e[1][0],i=+e[1][1]),n()):null==h?null:[[h,t],[r,i]]},scale:function(t){return arguments.length?(p=et((c=+t)*f,c*s,a,l),n()):c},translate:function(t){return arguments.length?(p=et(c*f,c*s,a=+t[0],l=+t[1]),n()):[a,l]},reflectX:function(t){
return arguments.length?(p=et(c*(f=t?-1:1),c*s,a,l),n()):f<0},reflectY:function(t){return arguments.length?(p=et(c*f,c*(s=t?-1:1),a,l),n()):s<0},fitExtent:function(n,t){return Fn(u,n,t)},fitSize:function(n,t){return In(u,n,t)}}};ot.invert=Qn(o);var Zi=function(){return Un(ot).scale(249.5).clipAngle(90+Bt)};ut.invert=Qn(function(n){return 2*Kt(n)});var Hi=function(){return Un(ut).scale(250).clipAngle(142)};ct.invert=function(n,t){return[-t,2*Kt($t(n))-Ut]};var Ji=function(){var n=Wn(ct),t=n.center,r=n.rotate;return n.center=function(n){return arguments.length?t([-n[1],n[0]]):(n=t(),[n[1],-n[0]])},n.rotate=function(n){return arguments.length?r([n[0],n[1],n.length>2?n[2]+90:90]):(n=r(),[n[0],n[1],n[2]-90])},r([0,0,90]).scale(159.155)};n.geoArea=pr,n.geoBounds=vr,n.geoCentroid=Er,n.geoCircle=zr,n.geoClipExtent=Fr,n.geoContains=Qr,n.geoDistance=Hr,n.geoGraticule=sn,n.geoGraticule10=pn,n.geoInterpolate=Vr,n.geoLength=Xr,n.geoPath=Ni,n.geoAlbers=Li,n.geoAlbersUsa=Oi,n.geoAzimuthalEqualArea=Ti,n.geoAzimuthalEqualAreaRaw=Gi,n.geoAzimuthalEquidistant=Fi,n.geoAzimuthalEquidistantRaw=ki,n.geoConicConformal=Bi,n.geoConicConformalRaw=nt,n.geoConicEqualArea=bi,n.geoConicEqualAreaRaw=Hn,n.geoConicEquidistant=Ui,n.geoConicEquidistantRaw=rt,n.geoEquirectangular=Di,n.geoEquirectangularRaw=tt,n.geoGnomonic=Xi,n.geoGnomonicRaw=it,n.geoIdentity=Yi,n.geoProjection=Un,n.geoProjectionMutator=Xn,n.geoMercator=Ii,n.geoMercatorRaw=Vn,n.geoOrthographic=Zi,n.geoOrthographicRaw=ot,n.geoStereographic=Hi,n.geoStereographicRaw=ut,n.geoTransverseMercator=Ji,n.geoTransverseMercatorRaw=ct,n.geoRotation=qr,n.geoStream=ar,n.geoTransform=Ai,Object.defineProperty(n,"__esModule",{value:!0})});');
INSERT INTO plv8_modules (modname, load_on_start, code) VALUES ('d3_contour', true, '// https://d3js.org/d3-contour/ Version 1.1.1. Copyright 2017 Mike Bostock.
!function(t,r){"object"==typeof exports&&"undefined"!=typeof module?r(exports,require("d3-array")):"function"==typeof define&&define.amd?define(["exports","d3-array"],r):r(t.d3=t.d3||{},t.d3)}(this,function(t,r){"use strict";function n(t,r){for(var n=r[0],i=r[1],a=-1,o=0,h=t.length,f=h-1;o<h;f=o++){var u=t[o],c=u[0],d=u[1],l=t[f],s=l[0],g=l[1];if(e(u,l,r))return 0;d>i!=g>i&&n<(s-c)*(i-d)/(g-d)+c&&(a=-a)}return a}function e(t,r,n){var e;return i(t,r,n)&&a(t[e=+(t[0]===r[0])],n[e],r[e])}function i(t,r,n){return(r[0]-t[0])*(n[1]-t[1])==(n[0]-t[0])*(r[1]-t[1])}function a(t,r,n){return t<=r&&r<=n||n<=r&&r<=t}function o(t,r,n){for(var e=t.width,i=t.height,a=1+(n<<1),o=0;o<i;++o)for(var h=0,f=0;h<e+n;++h)h<e&&(f+=t.data[h+o*e]),h>=n&&(h>=a&&(f-=t.data[h-a+o*e]),r.data[h-n+o*e]=f/Math.min(h+1,e-1+a-h,a))}function h(t,r,n){for(var e=t.width,i=t.height,a=1+(n<<1),o=0;o<e;++o)for(var h=0,f=0;h<i+n;++h)h<i&&(f+=t.data[o+h*e]),h>=n&&(h>=a&&(f-=t.data[o+(h-a)*e]),r.data[o+(h-n)*e]=f/Math.min(h+1,i-1+a-h,a))}function f(t){return t[0]}function u(t){return t[1]}var c=Array.prototype.slice,d=function(t,r){return t-r},l=function(t){for(var r=0,n=t.length,e=t[n-1][1]*t[0][0]-t[n-1][0]*t[0][1];++r<n;)e+=t[r-1][1]*t[r][0]-t[r-1][0]*t[r][1];return e},s=function(t){return function(){return t}},g=function(t,r){for(var e,i=-1,a=r.length;++i<a;)if(e=n(t,r[i]))return e;return 0},v=function(){},w=[[],[[[1,1.5],[.5,1]]],[[[1.5,1],[1,1.5]]],[[[1.5,1],[.5,1]]],[[[1,.5],[1.5,1]]],[[[1,.5],[.5,1]],[[1,1.5],[1.5,1]]],[[[1,.5],[1,1.5]]],[[[1,.5],[.5,1]]],[[[.5,1],[1,.5]]],[[[1,1.5],[1,.5]]],[[[.5,1],[1,1.5]],[[1.5,1],[1,.5]]],[[[1.5,1],[1,.5]]],[[[.5,1],[1.5,1]]],[[[1,1.5],[1.5,1]]],[[[.5,1],[1,1.5]]],[]],p=function(){function t(t){var e=h(t);if(Array.isArray(e))e=e.slice().sort(d);else{var i=r.extent(t),a=i[0],o=i[1];e=r.tickStep(a,o,e),e=r.range(Math.floor(a/e)*e,Math.floor(o/e)*e,e)}return e.map(function(r){var e=[],i=[];return n(t,r,function(n){f(n,t,r),l(n)>0?e.push([n]):i.push(n)}),i.forEach(function(t){for(var r,n=0,i=e.length;n<i;++n)if(-1!==g((r=e[n])[0],t))return void r.push(t)}),e}).map(function(t,r){return{type:"MultiPolygon",value:e[r],coordinates:t}})}function n(t,r,n){function i(t){var r,i,a=[t[0][0]+h,t[0][1]+f],o=[t[1][0]+h,t[1][1]+f],u=e(a),c=e(o);(r=g[u])?(i=s[c])?(delete g[r.end],delete s[i.start],r===i?(r.ring.push(o),n(r.ring)):s[r.start]=g[i.end]={start:r.start,end:i.end,ring:r.ring.concat(i.ring)}):(delete g[r.end],r.ring.push(o),g[r.end=c]=r):(r=s[c])?(i=g[u])?(delete s[r.start],delete g[i.end],r===i?(r.ring.push(o),n(r.ring)):s[i.start]=g[r.end]={start:i.start,end:r.end,ring:i.ring.concat(r.ring)}):(delete s[r.start],r.ring.unshift(a),s[r.start=u]=r):s[u]=g[c]={start:u,end:c,ring:[a,o]}}var h,f,u,c,d,l,s=new Array,g=new Array;for(h=f=-1,c=t[0]>=r,w[c<<1].forEach(i);++h<a-1;)u=c,c=t[h+1]>=r,w[u|c<<1].forEach(i);for(w[c<<0].forEach(i);++f<o-1;){for(h=-1,c=t[f*a+a]>=r,d=t[f*a]>=r,w[c<<1|d<<2].forEach(i);++h<a-1;)u=c,c=t[f*a+a+h+1]>=r,l=d,d=t[f*a+h+1]>=r,w[u|c<<1|d<<2|l<<3].forEach(i);w[c|d<<3].forEach(i)}for(h=-1,d=t[f*a]>=r,w[d<<2].forEach(i);++h<a-1;)l=d,d=t[f*a+h+1]>=r,w[d<<2|l<<3].forEach(i);w[d<<3].forEach(i)}function e(t){return 2*t[0]+t[1]*(a+1)*4}function i(t,r,n){t.forEach(function(t){var e,i=t[0],h=t[1],f=0|i,u=0|h,c=r[u*a+f];i>0&&i<a&&f===i&&(e=r[u*a+f-1],t[0]=i+(n-e)/(c-e)-.5),h>0&&h<o&&u===h&&(e=r[(u-1)*a+f],t[1]=h+(n-e)/(c-e)-.5)})}var a=1,o=1,h=r.thresholdSturges,f=i;return t.size=function(r){if(!arguments.length)return[a,o];var n=Math.ceil(r[0]),e=Math.ceil(r[1]);if(!(n>0&&e>0))throw new Error("invalid size");return a=n,o=e,t},t.thresholds=function(r){return arguments.length?(h="function"==typeof r?r:s(Array.isArray(r)?c.call(r):r),t):h},t.smooth=function(r){return arguments.length?(f=r?i:v,t):f===i},t};t.contours=p,t.contourDensity=function(){function t(t){var e=new Float32Array(A*m),i=new Float32Array(A*m);t.forEach(function(t,r,n){var i=l(t,r,n)+E>>M,a=g(t,r,n)+E>>M;i>=0&&i<A&&a>=0&&a<m&&++e[i+a*A]}),o({width:A,height:m,data:e},{width:A,height:m,data:i},y>>M),h({width:A,height:m,data:i},{width:A,height:m,data:e},y>>M),o({width:A,height:m,data:e},{width:A,height:m,data:i},y>>M),h({width:A,height:m,data:i},{width:A,height:m,data:e},y>>M),o({width:A,height:m,data:e},{width:A,height:m,data:i},y>>M),h({width:A,height:m,data:i},{width:A,height:m,data:e},y>>M);var a=z(e);if(!Array.isArray(a)){var f=r.max(e);a=r.tickStep(0,f,a),(a=r.range(0,Math.floor(f/a)*a,a)).shift()}return p().thresholds(a).size([A,m])(e).map(n)}function n(t){return t.value*=Math.pow(2,-2*M),t.coordinates.forEach(e),t}function e(t){t.forEach(i)}function i(t){t.forEach(a)}function a(t){t[0]=t[0]*Math.pow(2,M)-E,t[1]=t[1]*Math.pow(2,M)-E}function d(){return E=3*y,A=v+2*E>>M,m=w+2*E>>M,t}var l=f,g=u,v=960,w=500,y=20,M=2,E=3*y,A=v+2*E>>M,m=w+2*E>>M,z=s(20);return t.x=function(r){return arguments.length?(l="function"==typeof r?r:s(+r),t):l},t.y=function(r){return arguments.length?(g="function"==typeof r?r:s(+r),t):g},t.size=function(t){if(!arguments.length)return[v,w];var r=Math.ceil(t[0]),n=Math.ceil(t[1]);if(!(r>=0||r>=0))throw new Error("invalid size");return v=r,w=n,d()},t.cellSize=function(t){if(!arguments.length)return 1<<M;if(!((t=+t)>=1))throw new Error("invalid cell size");return M=Math.floor(Math.log(t)/Math.LN2),d()},t.thresholds=function(r){return arguments.length?(z="function"==typeof r?r:s(Array.isArray(r)?c.call(r):r),t):z},t.bandwidth=function(t){if(!arguments.length)return Math.sqrt(y*(y+1));if(!((t=+t)>=0))throw new Error("invalid bandwidth");return y=Math.round((Math.sqrt(4*t*t+1)-1)/2),d()},t},Object.defineProperty(t,"__esModule",{value:!0})});');
INSERT INTO plv8_modules (modname, load_on_start, code) VALUES ('d3_hexbin', true, '// https://github.com/d3/d3-hexbin Version 0.2.2. Copyright 2017 Mike Bostock.
!function(n,t){"object"==typeof exports&&"undefined"!=typeof module?t(exports):"function"==typeof define&&define.amd?define(["exports"],t):t(n.d3=n.d3||{})}(this,function(n){"use strict";function t(n){return n[0]}function r(n){return n[1]}var e=Math.PI/3,u=[0,e,2*e,3*e,4*e,5*e],o=function(){function n(n){var t,r={},e=[],u=n.length;for(t=0;t<u;++t)if(!isNaN(i=+d.call(null,o=n[t],t,n))&&!isNaN(c=+p.call(null,o,t,n))){var o,i,c,s=Math.round(c/=f),h=Math.round(i=i/a-(1&s)/2),l=c-s;if(3*Math.abs(l)>1){var v=i-h,M=h+(i<h?-1:1)/2,x=s+(c<s?-1:1),g=i-M,m=c-x;v*v+l*l>g*g+m*m&&(h=M+(1&s?1:-1)/2,s=x)}var y=h+"-"+s,j=r[y];j?j.push(o):(e.push(j=r[y]=[o]),j.x=(h+(1&s)/2)*a,j.y=s*f)}return e}function o(n){var t=0,r=0;return u.map(function(e){var u=Math.sin(e)*n,o=-Math.cos(e)*n,i=u-t,a=o-r;return t=u,r=o,[i,a]})}var i,a,f,c=0,s=0,h=1,l=1,d=t,p=r;return n.hexagon=function(n){return"m"+o(null==n?i:+n).join("l")+"z"},n.centers=function(){for(var n=[],t=Math.round(s/f),r=Math.round(c/a),e=t*f;e<l+i;e+=f,++t)for(var u=r*a+(1&t)*a/2;u<h+a/2;u+=a)n.push([u,e]);return n},n.mesh=function(){var t=o(i).slice(0,4).join("l");return n.centers().map(function(n){return"M"+n+"m"+t}).join("")},n.x=function(t){return arguments.length?(d=t,n):d},n.y=function(t){return arguments.length?(p=t,n):p},n.radius=function(t){return arguments.length?(i=+t,a=2*i*Math.sin(e),f=1.5*i,n):i},n.size=function(t){return arguments.length?(c=s=0,h=+t[0],l=+t[1],n):[h-c,l-s]},n.extent=function(t){return arguments.length?(c=+t[0][0],s=+t[0][1],h=+t[1][0],l=+t[1][1],n):[[c,s],[h,l]]},n.radius(1)};n.hexbin=o,Object.defineProperty(n,"__esModule",{value:!0})});');
INSERT INTO plv8_modules (modname, load_on_start, code) VALUES ('delaunator', true, '');
INSERT INTO plv8_modules (modname, load_on_start, code) VALUES ('topojson', true, '// https://github.com/topojson/topojson Version 3.0.0. Copyright 2017 Mike Bostock.
(function (global, factory) {
	typeof exports === ''object'' && typeof module !== ''undefined'' ? factory(exports) :
	typeof define === ''function'' && define.amd ? define([''exports''], factory) :
	(factory((global.topojson = global.topojson || {})));
}(this, (function (exports) { ''use strict'';

var identity = function(x) {
  return x;
};

var transform = function(transform) {
  if (transform == null) return identity;
  var x0,
      y0,
      kx = transform.scale[0],
      ky = transform.scale[1],
      dx = transform.translate[0],
      dy = transform.translate[1];
  return function(input, i) {
    if (!i) x0 = y0 = 0;
    var j = 2, n = input.length, output = new Array(n);
    output[0] = (x0 += input[0]) * kx + dx;
    output[1] = (y0 += input[1]) * ky + dy;
    while (j < n) output[j] = input[j], ++j;
    return output;
  };
};

var bbox = function(topology) {
  var t = transform(topology.transform), key,
      x0 = Infinity, y0 = x0, x1 = -x0, y1 = -x0;

  function bboxPoint(p) {
    p = t(p);
    if (p[0] < x0) x0 = p[0];
    if (p[0] > x1) x1 = p[0];
    if (p[1] < y0) y0 = p[1];
    if (p[1] > y1) y1 = p[1];
  }

  function bboxGeometry(o) {
    switch (o.type) {
      case "GeometryCollection": o.geometries.forEach(bboxGeometry); break;
      case "Point": bboxPoint(o.coordinates); break;
      case "MultiPoint": o.coordinates.forEach(bboxPoint); break;
    }
  }

  topology.arcs.forEach(function(arc) {
    var i = -1, n = arc.length, p;
    while (++i < n) {
      p = t(arc[i], i);
      if (p[0] < x0) x0 = p[0];
      if (p[0] > x1) x1 = p[0];
      if (p[1] < y0) y0 = p[1];
      if (p[1] > y1) y1 = p[1];
    }
  });

  for (key in topology.objects) {
    bboxGeometry(topology.objects[key]);
  }

  return [x0, y0, x1, y1];
};

var reverse = function(array, n) {
  var t, j = array.length, i = j - n;
  while (i < --j) t = array[i], array[i++] = array[j], array[j] = t;
};

var feature = function(topology, o) {
  return o.type === "GeometryCollection"
      ? {type: "FeatureCollection", features: o.geometries.map(function(o) { return feature$1(topology, o); })}
      : feature$1(topology, o);
};

function feature$1(topology, o) {
  var id = o.id,
      bbox = o.bbox,
      properties = o.properties == null ? {} : o.properties,
      geometry = object(topology, o);
  return id == null && bbox == null ? {type: "Feature", properties: properties, geometry: geometry}
      : bbox == null ? {type: "Feature", id: id, properties: properties, geometry: geometry}
      : {type: "Feature", id: id, bbox: bbox, properties: properties, geometry: geometry};
}

function object(topology, o) {
  var transformPoint = transform(topology.transform),
      arcs = topology.arcs;

  function arc(i, points) {
    if (points.length) points.pop();
    for (var a = arcs[i < 0 ? ~i : i], k = 0, n = a.length; k < n; ++k) {
      points.push(transformPoint(a[k], k));
    }
    if (i < 0) reverse(points, n);
  }

  function point(p) {
    return transformPoint(p);
  }

  function line(arcs) {
    var points = [];
    for (var i = 0, n = arcs.length; i < n; ++i) arc(arcs[i], points);
    if (points.length < 2) points.push(points[0]); // This should never happen per the specification.
    return points;
  }

  function ring(arcs) {
    var points = line(arcs);
    while (points.length < 4) points.push(points[0]); // This may happen if an arc has only two points.
    return points;
  }

  function polygon(arcs) {
    return arcs.map(ring);
  }

  function geometry(o) {
    var type = o.type, coordinates;
    switch (type) {
      case "GeometryCollection": return {type: type, geometries: o.geometries.map(geometry)};
      case "Point": coordinates = point(o.coordinates); break;
      case "MultiPoint": coordinates = o.coordinates.map(point); break;
      case "LineString": coordinates = line(o.arcs); break;
      case "MultiLineString": coordinates = o.arcs.map(line); break;
      case "Polygon": coordinates = polygon(o.arcs); break;
      case "MultiPolygon": coordinates = o.arcs.map(polygon); break;
      default: return null;
    }
    return {type: type, coordinates: coordinates};
  }

  return geometry(o);
}

var stitch = function(topology, arcs) {
  var stitchedArcs = {},
      fragmentByStart = {},
      fragmentByEnd = {},
      fragments = [],
      emptyIndex = -1;

  // Stitch empty arcs first, since they may be subsumed by other arcs.
  arcs.forEach(function(i, j) {
    var arc = topology.arcs[i < 0 ? ~i : i], t;
    if (arc.length < 3 && !arc[1][0] && !arc[1][1]) {
      t = arcs[++emptyIndex], arcs[emptyIndex] = i, arcs[j] = t;
    }
  });

  arcs.forEach(function(i) {
    var e = ends(i),
        start = e[0],
        end = e[1],
        f, g;

    if (f = fragmentByEnd[start]) {
      delete fragmentByEnd[f.end];
      f.push(i);
      f.end = end;
      if (g = fragmentByStart[end]) {
        delete fragmentByStart[g.start];
        var fg = g === f ? f : f.concat(g);
        fragmentByStart[fg.start = f.start] = fragmentByEnd[fg.end = g.end] = fg;
      } else {
        fragmentByStart[f.start] = fragmentByEnd[f.end] = f;
      }
    } else if (f = fragmentByStart[end]) {
      delete fragmentByStart[f.start];
      f.unshift(i);
      f.start = start;
      if (g = fragmentByEnd[start]) {
        delete fragmentByEnd[g.end];
        var gf = g === f ? f : g.concat(f);
        fragmentByStart[gf.start = g.start] = fragmentByEnd[gf.end = f.end] = gf;
      } else {
        fragmentByStart[f.start] = fragmentByEnd[f.end] = f;
      }
    } else {
      f = [i];
      fragmentByStart[f.start = start] = fragmentByEnd[f.end = end] = f;
    }
  });

  function ends(i) {
    var arc = topology.arcs[i < 0 ? ~i : i], p0 = arc[0], p1;
    if (topology.transform) p1 = [0, 0], arc.forEach(function(dp) { p1[0] += dp[0], p1[1] += dp[1]; });
    else p1 = arc[arc.length - 1];
    return i < 0 ? [p1, p0] : [p0, p1];
  }

  function flush(fragmentByEnd, fragmentByStart) {
    for (var k in fragmentByEnd) {
      var f = fragmentByEnd[k];
      delete fragmentByStart[f.start];
      delete f.start;
      delete f.end;
      f.forEach(function(i) { stitchedArcs[i < 0 ? ~i : i] = 1; });
      fragments.push(f);
    }
  }

  flush(fragmentByEnd, fragmentByStart);
  flush(fragmentByStart, fragmentByEnd);
  arcs.forEach(function(i) { if (!stitchedArcs[i < 0 ? ~i : i]) fragments.push([i]); });

  return fragments;
};

var mesh = function(topology) {
  return object(topology, meshArcs.apply(this, arguments));
};

function meshArcs(topology, object$$1, filter) {
  var arcs, i, n;
  if (arguments.length > 1) arcs = extractArcs(topology, object$$1, filter);
  else for (i = 0, arcs = new Array(n = topology.arcs.length); i < n; ++i) arcs[i] = i;
  return {type: "MultiLineString", arcs: stitch(topology, arcs)};
}

function extractArcs(topology, object$$1, filter) {
  var arcs = [],
      geomsByArc = [],
      geom;

  function extract0(i) {
    var j = i < 0 ? ~i : i;
    (geomsByArc[j] || (geomsByArc[j] = [])).push({i: i, g: geom});
  }

  function extract1(arcs) {
    arcs.forEach(extract0);
  }

  function extract2(arcs) {
    arcs.forEach(extract1);
  }

  function extract3(arcs) {
    arcs.forEach(extract2);
  }

  function geometry(o) {
    switch (geom = o, o.type) {
      case "GeometryCollection": o.geometries.forEach(geometry); break;
      case "LineString": extract1(o.arcs); break;
      case "MultiLineString": case "Polygon": extract2(o.arcs); break;
      case "MultiPolygon": extract3(o.arcs); break;
    }
  }

  geometry(object$$1);

  geomsByArc.forEach(filter == null
      ? function(geoms) { arcs.push(geoms[0].i); }
      : function(geoms) { if (filter(geoms[0].g, geoms[geoms.length - 1].g)) arcs.push(geoms[0].i); });

  return arcs;
}

function planarRingArea(ring) {
  var i = -1, n = ring.length, a, b = ring[n - 1], area = 0;
  while (++i < n) a = b, b = ring[i], area += a[0] * b[1] - a[1] * b[0];
  return Math.abs(area); // Note: doubled area!
}

var merge = function(topology) {
  return object(topology, mergeArcs.apply(this, arguments));
};

function mergeArcs(topology, objects) {
  var polygonsByArc = {},
      polygons = [],
      groups = [];

  objects.forEach(geometry);

  function geometry(o) {
    switch (o.type) {
      case "GeometryCollection": o.geometries.forEach(geometry); break;
      case "Polygon": extract(o.arcs); break;
      case "MultiPolygon": o.arcs.forEach(extract); break;
    }
  }

  function extract(polygon) {
    polygon.forEach(function(ring) {
      ring.forEach(function(arc) {
        (polygonsByArc[arc = arc < 0 ? ~arc : arc] || (polygonsByArc[arc] = [])).push(polygon);
      });
    });
    polygons.push(polygon);
  }

  function area(ring) {
    return planarRingArea(object(topology, {type: "Polygon", arcs: [ring]}).coordinates[0]);
  }

  polygons.forEach(function(polygon) {
    if (!polygon._) {
      var group = [],
          neighbors = [polygon];
      polygon._ = 1;
      groups.push(group);
      while (polygon = neighbors.pop()) {
        group.push(polygon);
        polygon.forEach(function(ring) {
          ring.forEach(function(arc) {
            polygonsByArc[arc < 0 ? ~arc : arc].forEach(function(polygon) {
              if (!polygon._) {
                polygon._ = 1;
                neighbors.push(polygon);
              }
            });
          });
        });
      }
    }
  });

  polygons.forEach(function(polygon) {
    delete polygon._;
  });

  return {
    type: "MultiPolygon",
    arcs: groups.map(function(polygons) {
      var arcs = [], n;

      // Extract the exterior (unique) arcs.
      polygons.forEach(function(polygon) {
        polygon.forEach(function(ring) {
          ring.forEach(function(arc) {
            if (polygonsByArc[arc < 0 ? ~arc : arc].length < 2) {
              arcs.push(arc);
            }
          });
        });
      });

      // Stitch the arcs into one or more rings.
      arcs = stitch(topology, arcs);

      // If more than one ring is returned,
      // at most one of these rings can be the exterior;
      // choose the one with the greatest absolute area.
      if ((n = arcs.length) > 1) {
        for (var i = 1, k = area(arcs[0]), ki, t; i < n; ++i) {
          if ((ki = area(arcs[i])) > k) {
            t = arcs[0], arcs[0] = arcs[i], arcs[i] = t, k = ki;
          }
        }
      }

      return arcs;
    })
  };
}

var bisect = function(a, x) {
  var lo = 0, hi = a.length;
  while (lo < hi) {
    var mid = lo + hi >>> 1;
    if (a[mid] < x) lo = mid + 1;
    else hi = mid;
  }
  return lo;
};

var neighbors = function(objects) {
  var indexesByArc = {}, // arc index -> array of object indexes
      neighbors = objects.map(function() { return []; });

  function line(arcs, i) {
    arcs.forEach(function(a) {
      if (a < 0) a = ~a;
      var o = indexesByArc[a];
      if (o) o.push(i);
      else indexesByArc[a] = [i];
    });
  }

  function polygon(arcs, i) {
    arcs.forEach(function(arc) { line(arc, i); });
  }

  function geometry(o, i) {
    if (o.type === "GeometryCollection") o.geometries.forEach(function(o) { geometry(o, i); });
    else if (o.type in geometryType) geometryType[o.type](o.arcs, i);
  }

  var geometryType = {
    LineString: line,
    MultiLineString: polygon,
    Polygon: polygon,
    MultiPolygon: function(arcs, i) { arcs.forEach(function(arc) { polygon(arc, i); }); }
  };

  objects.forEach(geometry);

  for (var i in indexesByArc) {
    for (var indexes = indexesByArc[i], m = indexes.length, j = 0; j < m; ++j) {
      for (var k = j + 1; k < m; ++k) {
        var ij = indexes[j], ik = indexes[k], n;
        if ((n = neighbors[ij])[i = bisect(n, ik)] !== ik) n.splice(i, 0, ik);
        if ((n = neighbors[ik])[i = bisect(n, ij)] !== ij) n.splice(i, 0, ij);
      }
    }
  }

  return neighbors;
};

var untransform = function(transform) {
  if (transform == null) return identity;
  var x0,
      y0,
      kx = transform.scale[0],
      ky = transform.scale[1],
      dx = transform.translate[0],
      dy = transform.translate[1];
  return function(input, i) {
    if (!i) x0 = y0 = 0;
    var j = 2,
        n = input.length,
        output = new Array(n),
        x1 = Math.round((input[0] - dx) / kx),
        y1 = Math.round((input[1] - dy) / ky);
    output[0] = x1 - x0, x0 = x1;
    output[1] = y1 - y0, y0 = y1;
    while (j < n) output[j] = input[j], ++j;
    return output;
  };
};

var quantize = function(topology, transform) {
  if (topology.transform) throw new Error("already quantized");

  if (!transform || !transform.scale) {
    if (!((n = Math.floor(transform)) >= 2)) throw new Error("n must be u22652");
    box = topology.bbox || bbox(topology);
    var x0 = box[0], y0 = box[1], x1 = box[2], y1 = box[3], n;
    transform = {scale: [x1 - x0 ? (x1 - x0) / (n - 1) : 1, y1 - y0 ? (y1 - y0) / (n - 1) : 1], translate: [x0, y0]};
  } else {
    box = topology.bbox;
  }

  var t = untransform(transform), box, key, inputs = topology.objects, outputs = {};

  function quantizePoint(point) {
    return t(point);
  }

  function quantizeGeometry(input) {
    var output;
    switch (input.type) {
      case "GeometryCollection": output = {type: "GeometryCollection", geometries: input.geometries.map(quantizeGeometry)}; break;
      case "Point": output = {type: "Point", coordinates: quantizePoint(input.coordinates)}; break;
      case "MultiPoint": output = {type: "MultiPoint", coordinates: input.coordinates.map(quantizePoint)}; break;
      default: return input;
    }
    if (input.id != null) output.id = input.id;
    if (input.bbox != null) output.bbox = input.bbox;
    if (input.properties != null) output.properties = input.properties;
    return output;
  }

  function quantizeArc(input) {
    var i = 0, j = 1, n = input.length, p, output = new Array(n); // pessimistic
    output[0] = t(input[0], 0);
    while (++i < n) if ((p = t(input[i], i))[0] || p[1]) output[j++] = p; // non-coincident points
    if (j === 1) output[j++] = [0, 0]; // an arc must have at least two points
    output.length = j;
    return output;
  }

  for (key in inputs) outputs[key] = quantizeGeometry(inputs[key]);

  return {
    type: "Topology",
    bbox: box,
    transform: transform,
    objects: outputs,
    arcs: topology.arcs.map(quantizeArc)
  };
};

// Computes the bounding box of the specified hash of GeoJSON objects.
var bounds = function(objects) {
  var x0 = Infinity,
      y0 = Infinity,
      x1 = -Infinity,
      y1 = -Infinity;

  function boundGeometry(geometry) {
    if (geometry != null && boundGeometryType.hasOwnProperty(geometry.type)) boundGeometryType[geometry.type](geometry);
  }

  var boundGeometryType = {
    GeometryCollection: function(o) { o.geometries.forEach(boundGeometry); },
    Point: function(o) { boundPoint(o.coordinates); },
    MultiPoint: function(o) { o.coordinates.forEach(boundPoint); },
    LineString: function(o) { boundLine(o.arcs); },
    MultiLineString: function(o) { o.arcs.forEach(boundLine); },
    Polygon: function(o) { o.arcs.forEach(boundLine); },
    MultiPolygon: function(o) { o.arcs.forEach(boundMultiLine); }
  };

  function boundPoint(coordinates) {
    var x = coordinates[0],
        y = coordinates[1];
    if (x < x0) x0 = x;
    if (x > x1) x1 = x;
    if (y < y0) y0 = y;
    if (y > y1) y1 = y;
  }

  function boundLine(coordinates) {
    coordinates.forEach(boundPoint);
  }

  function boundMultiLine(coordinates) {
    coordinates.forEach(boundLine);
  }

  for (var key in objects) {
    boundGeometry(objects[key]);
  }

  return x1 >= x0 && y1 >= y0 ? [x0, y0, x1, y1] : undefined;
};

var hashset = function(size, hash, equal, type, empty) {
  if (arguments.length === 3) {
    type = Array;
    empty = null;
  }

  var store = new type(size = 1 << Math.max(4, Math.ceil(Math.log(size) / Math.LN2))),
      mask = size - 1;

  for (var i = 0; i < size; ++i) {
    store[i] = empty;
  }

  function add(value) {
    var index = hash(value) & mask,
        match = store[index],
        collisions = 0;
    while (match != empty) {
      if (equal(match, value)) return true;
      if (++collisions >= size) throw new Error("full hashset");
      match = store[index = (index + 1) & mask];
    }
    store[index] = value;
    return true;
  }

  function has(value) {
    var index = hash(value) & mask,
        match = store[index],
        collisions = 0;
    while (match != empty) {
      if (equal(match, value)) return true;
      if (++collisions >= size) break;
      match = store[index = (index + 1) & mask];
    }
    return false;
  }

  function values() {
    var values = [];
    for (var i = 0, n = store.length; i < n; ++i) {
      var match = store[i];
      if (match != empty) values.push(match);
    }
    return values;
  }

  return {
    add: add,
    has: has,
    values: values
  };
};

var hashmap = function(size, hash, equal, keyType, keyEmpty, valueType) {
  if (arguments.length === 3) {
    keyType = valueType = Array;
    keyEmpty = null;
  }

  var keystore = new keyType(size = 1 << Math.max(4, Math.ceil(Math.log(size) / Math.LN2))),
      valstore = new valueType(size),
      mask = size - 1;

  for (var i = 0; i < size; ++i) {
    keystore[i] = keyEmpty;
  }

  function set(key, value) {
    var index = hash(key) & mask,
        matchKey = keystore[index],
        collisions = 0;
    while (matchKey != keyEmpty) {
      if (equal(matchKey, key)) return valstore[index] = value;
      if (++collisions >= size) throw new Error("full hashmap");
      matchKey = keystore[index = (index + 1) & mask];
    }
    keystore[index] = key;
    valstore[index] = value;
    return value;
  }

  function maybeSet(key, value) {
    var index = hash(key) & mask,
        matchKey = keystore[index],
        collisions = 0;
    while (matchKey != keyEmpty) {
      if (equal(matchKey, key)) return valstore[index];
      if (++collisions >= size) throw new Error("full hashmap");
      matchKey = keystore[index = (index + 1) & mask];
    }
    keystore[index] = key;
    valstore[index] = value;
    return value;
  }

  function get(key, missingValue) {
    var index = hash(key) & mask,
        matchKey = keystore[index],
        collisions = 0;
    while (matchKey != keyEmpty) {
      if (equal(matchKey, key)) return valstore[index];
      if (++collisions >= size) break;
      matchKey = keystore[index = (index + 1) & mask];
    }
    return missingValue;
  }

  function keys() {
    var keys = [];
    for (var i = 0, n = keystore.length; i < n; ++i) {
      var matchKey = keystore[i];
      if (matchKey != keyEmpty) keys.push(matchKey);
    }
    return keys;
  }

  return {
    set: set,
    maybeSet: maybeSet, // set if unset
    get: get,
    keys: keys
  };
};

var equalPoint = function(pointA, pointB) {
  return pointA[0] === pointB[0] && pointA[1] === pointB[1];
};

// TODO if quantized, use simpler Int32 hashing?

var buffer = new ArrayBuffer(16);
var floats = new Float64Array(buffer);
var uints = new Uint32Array(buffer);

var hashPoint = function(point) {
  floats[0] = point[0];
  floats[1] = point[1];
  var hash = uints[0] ^ uints[1];
  hash = hash << 5 ^ hash >> 7 ^ uints[2] ^ uints[3];
  return hash & 0x7fffffff;
};

// Given an extracted (pre-)topology, identifies all of the junctions. These are
// the points at which arcs (lines or rings) will need to be cut so that each
// arc is represented uniquely.
//
// A junction is a point where at least one arc deviates from another arc going
// through the same point. For example, consider the point B. If there is a arc
// through ABC and another arc through CBA, then B is not a junction because in
// both cases the adjacent point pairs are {A,C}. However, if there is an
// additional arc ABD, then {A,D} != {A,C}, and thus B becomes a junction.
//
// For a closed ring ABCA, the first point A’s adjacent points are the second
// and last point {B,C}. For a line, the first and last point are always
// considered junctions, even if the line is closed; this ensures that a closed
// line is never rotated.
var join = function(topology) {
  var coordinates = topology.coordinates,
      lines = topology.lines,
      rings = topology.rings,
      indexes = index(),
      visitedByIndex = new Int32Array(coordinates.length),
      leftByIndex = new Int32Array(coordinates.length),
      rightByIndex = new Int32Array(coordinates.length),
      junctionByIndex = new Int8Array(coordinates.length),
      junctionCount = 0, // upper bound on number of junctions
      i, n,
      previousIndex,
      currentIndex,
      nextIndex;

  for (i = 0, n = coordinates.length; i < n; ++i) {
    visitedByIndex[i] = leftByIndex[i] = rightByIndex[i] = -1;
  }

  for (i = 0, n = lines.length; i < n; ++i) {
    var line = lines[i],
        lineStart = line[0],
        lineEnd = line[1];
    currentIndex = indexes[lineStart];
    nextIndex = indexes[++lineStart];
    ++junctionCount, junctionByIndex[currentIndex] = 1; // start
    while (++lineStart <= lineEnd) {
      sequence(i, previousIndex = currentIndex, currentIndex = nextIndex, nextIndex = indexes[lineStart]);
    }
    ++junctionCount, junctionByIndex[nextIndex] = 1; // end
  }

  for (i = 0, n = coordinates.length; i < n; ++i) {
    visitedByIndex[i] = -1;
  }

  for (i = 0, n = rings.length; i < n; ++i) {
    var ring = rings[i],
        ringStart = ring[0] + 1,
        ringEnd = ring[1];
    previousIndex = indexes[ringEnd - 1];
    currentIndex = indexes[ringStart - 1];
    nextIndex = indexes[ringStart];
    sequence(i, previousIndex, currentIndex, nextIndex);
    while (++ringStart <= ringEnd) {
      sequence(i, previousIndex = currentIndex, currentIndex = nextIndex, nextIndex = indexes[ringStart]);
    }
  }

  function sequence(i, previousIndex, currentIndex, nextIndex) {
    if (visitedByIndex[currentIndex] === i) return; // ignore self-intersection
    visitedByIndex[currentIndex] = i;
    var leftIndex = leftByIndex[currentIndex];
    if (leftIndex >= 0) {
      var rightIndex = rightByIndex[currentIndex];
      if ((leftIndex !== previousIndex || rightIndex !== nextIndex)
        && (leftIndex !== nextIndex || rightIndex !== previousIndex)) {
        ++junctionCount, junctionByIndex[currentIndex] = 1;
      }
    } else {
      leftByIndex[currentIndex] = previousIndex;
      rightByIndex[currentIndex] = nextIndex;
    }
  }

  function index() {
    var indexByPoint = hashmap(coordinates.length * 1.4, hashIndex, equalIndex, Int32Array, -1, Int32Array),
        indexes = new Int32Array(coordinates.length);

    for (var i = 0, n = coordinates.length; i < n; ++i) {
      indexes[i] = indexByPoint.maybeSet(i, i);
    }

    return indexes;
  }

  function hashIndex(i) {
    return hashPoint(coordinates[i]);
  }

  function equalIndex(i, j) {
    return equalPoint(coordinates[i], coordinates[j]);
  }

  visitedByIndex = leftByIndex = rightByIndex = null;

  var junctionByPoint = hashset(junctionCount * 1.4, hashPoint, equalPoint), j;

  // Convert back to a standard hashset by point for caller convenience.
  for (i = 0, n = coordinates.length; i < n; ++i) {
    if (junctionByIndex[j = indexes[i]]) {
      junctionByPoint.add(coordinates[j]);
    }
  }

  return junctionByPoint;
};

// Given an extracted (pre-)topology, cuts (or rotates) arcs so that all shared
// point sequences are identified. The topology can then be subsequently deduped
// to remove exact duplicate arcs.
var cut = function(topology) {
  var junctions = join(topology),
      coordinates = topology.coordinates,
      lines = topology.lines,
      rings = topology.rings,
      next,
      i, n;

  for (i = 0, n = lines.length; i < n; ++i) {
    var line = lines[i],
        lineMid = line[0],
        lineEnd = line[1];
    while (++lineMid < lineEnd) {
      if (junctions.has(coordinates[lineMid])) {
        next = {0: lineMid, 1: line[1]};
        line[1] = lineMid;
        line = line.next = next;
      }
    }
  }

  for (i = 0, n = rings.length; i < n; ++i) {
    var ring = rings[i],
        ringStart = ring[0],
        ringMid = ringStart,
        ringEnd = ring[1],
        ringFixed = junctions.has(coordinates[ringStart]);
    while (++ringMid < ringEnd) {
      if (junctions.has(coordinates[ringMid])) {
        if (ringFixed) {
          next = {0: ringMid, 1: ring[1]};
          ring[1] = ringMid;
          ring = ring.next = next;
        } else { // For the first junction, we can rotate rather than cut.
          rotateArray(coordinates, ringStart, ringEnd, ringEnd - ringMid);
          coordinates[ringEnd] = coordinates[ringStart];
          ringFixed = true;
          ringMid = ringStart; // restart; we may have skipped junctions
        }
      }
    }
  }

  return topology;
};

function rotateArray(array, start, end, offset) {
  reverse$1(array, start, end);
  reverse$1(array, start, start + offset);
  reverse$1(array, start + offset, end);
}

function reverse$1(array, start, end) {
  for (var mid = start + ((end-- - start) >> 1), t; start < mid; ++start, --end) {
    t = array[start], array[start] = array[end], array[end] = t;
  }
}

// Given a cut topology, combines duplicate arcs.
var dedup = function(topology) {
  var coordinates = topology.coordinates,
      lines = topology.lines, line,
      rings = topology.rings, ring,
      arcCount = lines.length + rings.length,
      i, n;

  delete topology.lines;
  delete topology.rings;

  // Count the number of (non-unique) arcs to initialize the hashmap safely.
  for (i = 0, n = lines.length; i < n; ++i) {
    line = lines[i]; while (line = line.next) ++arcCount;
  }
  for (i = 0, n = rings.length; i < n; ++i) {
    ring = rings[i]; while (ring = ring.next) ++arcCount;
  }

  var arcsByEnd = hashmap(arcCount * 2 * 1.4, hashPoint, equalPoint),
      arcs = topology.arcs = [];

  for (i = 0, n = lines.length; i < n; ++i) {
    line = lines[i];
    do {
      dedupLine(line);
    } while (line = line.next);
  }

  for (i = 0, n = rings.length; i < n; ++i) {
    ring = rings[i];
    if (ring.next) { // arc is no longer closed
      do {
        dedupLine(ring);
      } while (ring = ring.next);
    } else {
      dedupRing(ring);
    }
  }

  function dedupLine(arc) {
    var startPoint,
        endPoint,
        startArcs, startArc,
        endArcs, endArc,
        i, n;

    // Does this arc match an existing arc in order?
    if (startArcs = arcsByEnd.get(startPoint = coordinates[arc[0]])) {
      for (i = 0, n = startArcs.length; i < n; ++i) {
        startArc = startArcs[i];
        if (equalLine(startArc, arc)) {
          arc[0] = startArc[0];
          arc[1] = startArc[1];
          return;
        }
      }
    }

    // Does this arc match an existing arc in reverse order?
    if (endArcs = arcsByEnd.get(endPoint = coordinates[arc[1]])) {
      for (i = 0, n = endArcs.length; i < n; ++i) {
        endArc = endArcs[i];
        if (reverseEqualLine(endArc, arc)) {
          arc[1] = endArc[0];
          arc[0] = endArc[1];
          return;
        }
      }
    }

    if (startArcs) startArcs.push(arc); else arcsByEnd.set(startPoint, [arc]);
    if (endArcs) endArcs.push(arc); else arcsByEnd.set(endPoint, [arc]);
    arcs.push(arc);
  }

  function dedupRing(arc) {
    var endPoint,
        endArcs,
        endArc,
        i, n;

    // Does this arc match an existing line in order, or reverse order?
    // Rings are closed, so their start point and end point is the same.
    if (endArcs = arcsByEnd.get(endPoint = coordinates[arc[0]])) {
      for (i = 0, n = endArcs.length; i < n; ++i) {
        endArc = endArcs[i];
        if (equalRing(endArc, arc)) {
          arc[0] = endArc[0];
          arc[1] = endArc[1];
          return;
        }
        if (reverseEqualRing(endArc, arc)) {
          arc[0] = endArc[1];
          arc[1] = endArc[0];
          return;
        }
      }
    }

    // Otherwise, does this arc match an existing ring in order, or reverse order?
    if (endArcs = arcsByEnd.get(endPoint = coordinates[arc[0] + findMinimumOffset(arc)])) {
      for (i = 0, n = endArcs.length; i < n; ++i) {
        endArc = endArcs[i];
        if (equalRing(endArc, arc)) {
          arc[0] = endArc[0];
          arc[1] = endArc[1];
          return;
        }
        if (reverseEqualRing(endArc, arc)) {
          arc[0] = endArc[1];
          arc[1] = endArc[0];
          return;
        }
      }
    }

    if (endArcs) endArcs.push(arc); else arcsByEnd.set(endPoint, [arc]);
    arcs.push(arc);
  }

  function equalLine(arcA, arcB) {
    var ia = arcA[0], ib = arcB[0],
        ja = arcA[1], jb = arcB[1];
    if (ia - ja !== ib - jb) return false;
    for (; ia <= ja; ++ia, ++ib) if (!equalPoint(coordinates[ia], coordinates[ib])) return false;
    return true;
  }

  function reverseEqualLine(arcA, arcB) {
    var ia = arcA[0], ib = arcB[0],
        ja = arcA[1], jb = arcB[1];
    if (ia - ja !== ib - jb) return false;
    for (; ia <= ja; ++ia, --jb) if (!equalPoint(coordinates[ia], coordinates[jb])) return false;
    return true;
  }

  function equalRing(arcA, arcB) {
    var ia = arcA[0], ib = arcB[0],
        ja = arcA[1], jb = arcB[1],
        n = ja - ia;
    if (n !== jb - ib) return false;
    var ka = findMinimumOffset(arcA),
        kb = findMinimumOffset(arcB);
    for (var i = 0; i < n; ++i) {
      if (!equalPoint(coordinates[ia + (i + ka) % n], coordinates[ib + (i + kb) % n])) return false;
    }
    return true;
  }

  function reverseEqualRing(arcA, arcB) {
    var ia = arcA[0], ib = arcB[0],
        ja = arcA[1], jb = arcB[1],
        n = ja - ia;
    if (n !== jb - ib) return false;
    var ka = findMinimumOffset(arcA),
        kb = n - findMinimumOffset(arcB);
    for (var i = 0; i < n; ++i) {
      if (!equalPoint(coordinates[ia + (i + ka) % n], coordinates[jb - (i + kb) % n])) return false;
    }
    return true;
  }

  // Rings are rotated to a consistent, but arbitrary, start point.
  // This is necessary to detect when a ring and a rotated copy are dupes.
  function findMinimumOffset(arc) {
    var start = arc[0],
        end = arc[1],
        mid = start,
        minimum = mid,
        minimumPoint = coordinates[mid];
    while (++mid < end) {
      var point = coordinates[mid];
      if (point[0] < minimumPoint[0] || point[0] === minimumPoint[0] && point[1] < minimumPoint[1]) {
        minimum = mid;
        minimumPoint = point;
      }
    }
    return minimum - start;
  }

  return topology;
};

// Given an array of arcs in absolute (but already quantized!) coordinates,
// converts to fixed-point delta encoding.
// This is a destructive operation that modifies the given arcs!
var delta = function(arcs) {
  var i = -1,
      n = arcs.length;

  while (++i < n) {
    var arc = arcs[i],
        j = 0,
        k = 1,
        m = arc.length,
        point = arc[0],
        x0 = point[0],
        y0 = point[1],
        x1,
        y1;

    while (++j < m) {
      point = arc[j], x1 = point[0], y1 = point[1];
      if (x1 !== x0 || y1 !== y0) arc[k++] = [x1 - x0, y1 - y0], x0 = x1, y0 = y1;
    }

    if (k === 1) arc[k++] = [0, 0]; // Each arc must be an array of two or more positions.

    arc.length = k;
  }

  return arcs;
};

// Extracts the lines and rings from the specified hash of geometry objects.
//
// Returns an object with three properties:
//
// * coordinates - shared buffer of [x, y] coordinates
// * lines - lines extracted from the hash, of the form [start, end]
// * rings - rings extracted from the hash, of the form [start, end]
//
// For each ring or line, start and end represent inclusive indexes into the
// coordinates buffer. For rings (and closed lines), coordinates[start] equals
// coordinates[end].
//
// For each line or polygon geometry in the input hash, including nested
// geometries as in geometry collections, the `coordinates` array is replaced
// with an equivalent `arcs` array that, for each line (for line string
// geometries) or ring (for polygon geometries), points to one of the above
// lines or rings.
var extract = function(objects) {
  var index = -1,
      lines = [],
      rings = [],
      coordinates = [];

  function extractGeometry(geometry) {
    if (geometry && extractGeometryType.hasOwnProperty(geometry.type)) extractGeometryType[geometry.type](geometry);
  }

  var extractGeometryType = {
    GeometryCollection: function(o) { o.geometries.forEach(extractGeometry); },
    LineString: function(o) { o.arcs = extractLine(o.arcs); },
    MultiLineString: function(o) { o.arcs = o.arcs.map(extractLine); },
    Polygon: function(o) { o.arcs = o.arcs.map(extractRing); },
    MultiPolygon: function(o) { o.arcs = o.arcs.map(extractMultiRing); }
  };

  function extractLine(line) {
    for (var i = 0, n = line.length; i < n; ++i) coordinates[++index] = line[i];
    var arc = {0: index - n + 1, 1: index};
    lines.push(arc);
    return arc;
  }

  function extractRing(ring) {
    for (var i = 0, n = ring.length; i < n; ++i) coordinates[++index] = ring[i];
    var arc = {0: index - n + 1, 1: index};
    rings.push(arc);
    return arc;
  }

  function extractMultiRing(rings) {
    return rings.map(extractRing);
  }

  for (var key in objects) {
    extractGeometry(objects[key]);
  }

  return {
    type: "Topology",
    coordinates: coordinates,
    lines: lines,
    rings: rings,
    objects: objects
  };
};

// Given a hash of GeoJSON objects, returns a hash of GeoJSON geometry objects.
// Any null input geometry objects are represented as {type: null} in the output.
// Any feature.{id,properties,bbox} are transferred to the output geometry object.
// Each output geometry object is a shallow copy of the input (e.g., properties, coordinates)!
var geometry = function(inputs) {
  var outputs = {}, key;
  for (key in inputs) outputs[key] = geomifyObject(inputs[key]);
  return outputs;
};

function geomifyObject(input) {
  return input == null ? {type: null}
      : (input.type === "FeatureCollection" ? geomifyFeatureCollection
      : input.type === "Feature" ? geomifyFeature
      : geomifyGeometry)(input);
}

function geomifyFeatureCollection(input) {
  var output = {type: "GeometryCollection", geometries: input.features.map(geomifyFeature)};
  if (input.bbox != null) output.bbox = input.bbox;
  return output;
}

function geomifyFeature(input) {
  var output = geomifyGeometry(input.geometry), key; // eslint-disable-line no-unused-vars
  if (input.id != null) output.id = input.id;
  if (input.bbox != null) output.bbox = input.bbox;
  for (key in input.properties) { output.properties = input.properties; break; }
  return output;
}

function geomifyGeometry(input) {
  if (input == null) return {type: null};
  var output = input.type === "GeometryCollection" ? {type: "GeometryCollection", geometries: input.geometries.map(geomifyGeometry)}
      : input.type === "Point" || input.type === "MultiPoint" ? {type: input.type, coordinates: input.coordinates}
      : {type: input.type, arcs: input.coordinates}; // TODO Check for unknown types?
  if (input.bbox != null) output.bbox = input.bbox;
  return output;
}

var prequantize = function(objects, bbox, n) {
  var x0 = bbox[0],
      y0 = bbox[1],
      x1 = bbox[2],
      y1 = bbox[3],
      kx = x1 - x0 ? (n - 1) / (x1 - x0) : 1,
      ky = y1 - y0 ? (n - 1) / (y1 - y0) : 1;

  function quantizePoint(input) {
    return [Math.round((input[0] - x0) * kx), Math.round((input[1] - y0) * ky)];
  }

  function quantizePoints(input, m) {
    var i = -1,
        j = 0,
        n = input.length,
        output = new Array(n), // pessimistic
        pi,
        px,
        py,
        x,
        y;

    while (++i < n) {
      pi = input[i];
      x = Math.round((pi[0] - x0) * kx);
      y = Math.round((pi[1] - y0) * ky);
      if (x !== px || y !== py) output[j++] = [px = x, py = y]; // non-coincident points
    }

    output.length = j;
    while (j < m) j = output.push([output[0][0], output[0][1]]);
    return output;
  }

  function quantizeLine(input) {
    return quantizePoints(input, 2);
  }

  function quantizeRing(input) {
    return quantizePoints(input, 4);
  }

  function quantizePolygon(input) {
    return input.map(quantizeRing);
  }

  function quantizeGeometry(o) {
    if (o != null && quantizeGeometryType.hasOwnProperty(o.type)) quantizeGeometryType[o.type](o);
  }

  var quantizeGeometryType = {
    GeometryCollection: function(o) { o.geometries.forEach(quantizeGeometry); },
    Point: function(o) { o.coordinates = quantizePoint(o.coordinates); },
    MultiPoint: function(o) { o.coordinates = o.coordinates.map(quantizePoint); },
    LineString: function(o) { o.arcs = quantizeLine(o.arcs); },
    MultiLineString: function(o) { o.arcs = o.arcs.map(quantizeLine); },
    Polygon: function(o) { o.arcs = quantizePolygon(o.arcs); },
    MultiPolygon: function(o) { o.arcs = o.arcs.map(quantizePolygon); }
  };

  for (var key in objects) {
    quantizeGeometry(objects[key]);
  }

  return {
    scale: [1 / kx, 1 / ky],
    translate: [x0, y0]
  };
};

// Constructs the TopoJSON Topology for the specified hash of features.
// Each object in the specified hash must be a GeoJSON object,
// meaning FeatureCollection, a Feature or a geometry object.
var topology = function(objects, quantization) {
  var bbox = bounds(objects = geometry(objects)),
      transform = quantization > 0 && bbox && prequantize(objects, bbox, quantization),
      topology = dedup(cut(extract(objects))),
      coordinates = topology.coordinates,
      indexByArc = hashmap(topology.arcs.length * 1.4, hashArc, equalArc);

  objects = topology.objects; // for garbage collection
  topology.bbox = bbox;
  topology.arcs = topology.arcs.map(function(arc, i) {
    indexByArc.set(arc, i);
    return coordinates.slice(arc[0], arc[1] + 1);
  });

  delete topology.coordinates;
  coordinates = null;

  function indexGeometry(geometry$$1) {
    if (geometry$$1 && indexGeometryType.hasOwnProperty(geometry$$1.type)) indexGeometryType[geometry$$1.type](geometry$$1);
  }

  var indexGeometryType = {
    GeometryCollection: function(o) { o.geometries.forEach(indexGeometry); },
    LineString: function(o) { o.arcs = indexArcs(o.arcs); },
    MultiLineString: function(o) { o.arcs = o.arcs.map(indexArcs); },
    Polygon: function(o) { o.arcs = o.arcs.map(indexArcs); },
    MultiPolygon: function(o) { o.arcs = o.arcs.map(indexMultiArcs); }
  };

  function indexArcs(arc) {
    var indexes = [];
    do {
      var index = indexByArc.get(arc);
      indexes.push(arc[0] < arc[1] ? index : ~index);
    } while (arc = arc.next);
    return indexes;
  }

  function indexMultiArcs(arcs) {
    return arcs.map(indexArcs);
  }

  for (var key in objects) {
    indexGeometry(objects[key]);
  }

  if (transform) {
    topology.transform = transform;
    topology.arcs = delta(topology.arcs);
  }

  return topology;
};

function hashArc(arc) {
  var i = arc[0], j = arc[1], t;
  if (j < i) t = i, i = j, j = t;
  return i + 31 * j;
}

function equalArc(arcA, arcB) {
  var ia = arcA[0], ja = arcA[1],
      ib = arcB[0], jb = arcB[1], t;
  if (ja < ia) t = ia, ia = ja, ja = t;
  if (jb < ib) t = ib, ib = jb, jb = t;
  return ia === ib && ja === jb;
}

var prune = function(topology) {
  var oldObjects = topology.objects,
      newObjects = {},
      oldArcs = topology.arcs,
      newArcs = [],
      newArcIndex = -1,
      newIndexByOldIndex = new Array(oldArcs.length),
      key;

  function pruneGeometry(input) {
    var output;
    switch (input.type) {
      case "GeometryCollection": output = {type: "GeometryCollection", geometries: input.geometries.map(pruneGeometry)}; break;
      case "LineString": output = {type: "LineString", arcs: pruneArcs(input.arcs)}; break;
      case "MultiLineString": output = {type: "MultiLineString", arcs: input.arcs.map(pruneArcs)}; break;
      case "Polygon": output = {type: "Polygon", arcs: input.arcs.map(pruneArcs)}; break;
      case "MultiPolygon": output = {type: "MultiPolygon", arcs: input.arcs.map(pruneMultiArcs)}; break;
      default: return input;
    }
    if (input.id != null) output.id = input.id;
    if (input.bbox != null) output.bbox = input.bbox;
    if (input.properties != null) output.properties = input.properties;
    return output;
  }

  function pruneArc(oldIndex) {
    var oldReverse = oldIndex < 0 && (oldIndex = ~oldIndex, true), newIndex;

    // If this is the first instance of this arc, record it under its new index.
    if ((newIndex = newIndexByOldIndex[oldIndex]) == null) {
      newIndexByOldIndex[oldIndex] = newIndex = ++newArcIndex;
      newArcs[newIndex] = oldArcs[oldIndex];
    }

    return oldReverse ? ~newIndex : newIndex;
  }

  function pruneArcs(arcs) {
    return arcs.map(pruneArc);
  }

  function pruneMultiArcs(arcs) {
    return arcs.map(pruneArcs);
  }

  for (key in oldObjects) {
    newObjects[key] = pruneGeometry(oldObjects[key]);
  }

  return {
    type: "Topology",
    bbox: topology.bbox,
    transform: topology.transform,
    objects: newObjects,
    arcs: newArcs
  };
};

var filter = function(topology, filter) {
  var oldObjects = topology.objects,
      newObjects = {},
      key;

  if (filter == null) filter = filterTrue;

  function filterGeometry(input) {
    var output, arcs;
    switch (input.type) {
      case "Polygon": {
        arcs = filterRings(input.arcs);
        output = arcs ? {type: "Polygon", arcs: arcs} : {type: null};
        break;
      }
      case "MultiPolygon": {
        arcs = input.arcs.map(filterRings).filter(filterIdentity);
        output = arcs.length ? {type: "MultiPolygon", arcs: arcs} : {type: null};
        break;
      }
      case "GeometryCollection": {
        arcs = input.geometries.map(filterGeometry).filter(filterNotNull);
        output = arcs.length ? {type: "GeometryCollection", geometries: arcs} : {type: null};
        break;
      }
      default: return input;
    }
    if (input.id != null) output.id = input.id;
    if (input.bbox != null) output.bbox = input.bbox;
    if (input.properties != null) output.properties = input.properties;
    return output;
  }

  function filterRings(arcs) {
    return arcs.length && filterExteriorRing(arcs[0]) // if the exterior is small, ignore any holes
        ? [arcs[0]].concat(arcs.slice(1).filter(filterInteriorRing))
        : null;
  }

  function filterExteriorRing(ring) {
    return filter(ring, false);
  }

  function filterInteriorRing(ring) {
    return filter(ring, true);
  }

  for (key in oldObjects) {
    newObjects[key] = filterGeometry(oldObjects[key]);
  }

  return prune({
    type: "Topology",
    bbox: topology.bbox,
    transform: topology.transform,
    objects: newObjects,
    arcs: topology.arcs
  });
};

function filterTrue() {
  return true;
}

function filterIdentity(x) {
  return x;
}

function filterNotNull(geometry) {
  return geometry.type != null;
}

var filterAttached = function(topology) {
  var uniqueRingByArc = {}, // arc index -> index of unique associated ring, or -1 if used by multiple rings
      ringIndex = 0,
      name;

  function testGeometry(o) {
    switch (o.type) {
      case "GeometryCollection": o.geometries.forEach(testGeometry); break;
      case "Polygon": testArcs(o.arcs); break;
      case "MultiPolygon": o.arcs.forEach(testArcs); break;
    }
  }

  function testArcs(arcs) {
    for (var i = 0, n = arcs.length; i < n; ++i, ++ringIndex) {
      for (var ring = arcs[i], j = 0, m = ring.length; j < m; ++j) {
        var arc = ring[j];
        if (arc < 0) arc = ~arc;
        var uniqueRing = uniqueRingByArc[arc];
        if (uniqueRing >= 0 && uniqueRing !== ringIndex) uniqueRingByArc[arc] = -1;
        else uniqueRingByArc[arc] = ringIndex;
      }
    }
  }

  for (name in topology.objects) {
    testGeometry(topology.objects[name]);
  }

  return function(ring) {
    for (var j = 0, m = ring.length, arc; j < m; ++j) {
      if (arc = ring[j], uniqueRingByArc[arc < 0 ? ~arc : arc] < 0) {
        return true;
      }
    }
    return false;
  };
};

function planarTriangleArea(triangle) {
  var a = triangle[0], b = triangle[1], c = triangle[2];
  return Math.abs((a[0] - c[0]) * (b[1] - a[1]) - (a[0] - b[0]) * (c[1] - a[1])) / 2;
}

function planarRingArea$1(ring) {
  var i = -1, n = ring.length, a, b = ring[n - 1], area = 0;
  while (++i < n) a = b, b = ring[i], area += a[0] * b[1] - a[1] * b[0];
  return Math.abs(area) / 2;
}

var filterWeight = function(topology, minWeight, weight) {
  minWeight = minWeight == null ? Number.MIN_VALUE : +minWeight;

  if (weight == null) weight = planarRingArea$1;

  return function(ring, interior) {
    return weight(feature(topology, {type: "Polygon", arcs: [ring]}).geometry.coordinates[0], interior) >= minWeight;
  };
};

var filterAttachedWeight = function(topology, minWeight, weight) {
  var a = filterAttached(topology),
      w = filterWeight(topology, minWeight, weight);
  return function(ring, interior) {
    return a(ring, interior) || w(ring, interior);
  };
};

function compare(a, b) {
  return a[1][2] - b[1][2];
}

var newHeap = function() {
  var heap = {},
      array = [],
      size = 0;

  heap.push = function(object) {
    up(array[object._ = size] = object, size++);
    return size;
  };

  heap.pop = function() {
    if (size <= 0) return;
    var removed = array[0], object;
    if (--size > 0) object = array[size], down(array[object._ = 0] = object, 0);
    return removed;
  };

  heap.remove = function(removed) {
    var i = removed._, object;
    if (array[i] !== removed) return; // invalid request
    if (i !== --size) object = array[size], (compare(object, removed) < 0 ? up : down)(array[object._ = i] = object, i);
    return i;
  };

  function up(object, i) {
    while (i > 0) {
      var j = ((i + 1) >> 1) - 1,
          parent = array[j];
      if (compare(object, parent) >= 0) break;
      array[parent._ = i] = parent;
      array[object._ = i = j] = object;
    }
  }

  function down(object, i) {
    while (true) {
      var r = (i + 1) << 1,
          l = r - 1,
          j = i,
          child = array[j];
      if (l < size && compare(array[l], child) < 0) child = array[j = l];
      if (r < size && compare(array[r], child) < 0) child = array[j = r];
      if (j === i) break;
      array[child._ = i] = child;
      array[object._ = i = j] = object;
    }
  }

  return heap;
};

function copy(point) {
  return [point[0], point[1], 0];
}

var presimplify = function(topology, weight) {
  var point = topology.transform ? transform(topology.transform) : copy,
      heap = newHeap();

  if (weight == null) weight = planarTriangleArea;

  var arcs = topology.arcs.map(function(arc) {
    var triangles = [],
        maxWeight = 0,
        triangle,
        i,
        n;

    arc = arc.map(point);

    for (i = 1, n = arc.length - 1; i < n; ++i) {
      triangle = [arc[i - 1], arc[i], arc[i + 1]];
      triangle[1][2] = weight(triangle);
      triangles.push(triangle);
      heap.push(triangle);
    }

    // Always keep the arc endpoints!
    arc[0][2] = arc[n][2] = Infinity;

    for (i = 0, n = triangles.length; i < n; ++i) {
      triangle = triangles[i];
      triangle.previous = triangles[i - 1];
      triangle.next = triangles[i + 1];
    }

    while (triangle = heap.pop()) {
      var previous = triangle.previous,
          next = triangle.next;

      // If the weight of the current point is less than that of the previous
      // point to be eliminated, use the latter’s weight instead. This ensures
      // that the current point cannot be eliminated without eliminating
      // previously- eliminated points.
      if (triangle[1][2] < maxWeight) triangle[1][2] = maxWeight;
      else maxWeight = triangle[1][2];

      if (previous) {
        previous.next = next;
        previous[2] = triangle[2];
        update(previous);
      }

      if (next) {
        next.previous = previous;
        next[0] = triangle[0];
        update(next);
      }
    }

    return arc;
  });

  function update(triangle) {
    heap.remove(triangle);
    triangle[1][2] = weight(triangle);
    heap.push(triangle);
  }

  return {
    type: "Topology",
    bbox: topology.bbox,
    objects: topology.objects,
    arcs: arcs
  };
};

var quantile = function(topology, p) {
  var array = [];

  topology.arcs.forEach(function(arc) {
    arc.forEach(function(point) {
      if (isFinite(point[2])) { // Ignore endpoints, whose weight is Infinity.
        array.push(point[2]);
      }
    });
  });

  return array.length && quantile$1(array.sort(descending), p);
};

function quantile$1(array, p) {
  if (!(n = array.length)) return;
  if ((p = +p) <= 0 || n < 2) return array[0];
  if (p >= 1) return array[n - 1];
  var n,
      h = (n - 1) * p,
      i = Math.floor(h),
      a = array[i],
      b = array[i + 1];
  return a + (b - a) * (h - i);
}

function descending(a, b) {
  return b - a;
}

var simplify = function(topology, minWeight) {
  minWeight = minWeight == null ? Number.MIN_VALUE : +minWeight;

  // Remove points whose weight is less than the minimum weight.
  var arcs = topology.arcs.map(function(input) {
    var i = -1,
        j = 0,
        n = input.length,
        output = new Array(n), // pessimistic
        point;

    while (++i < n) {
      if ((point = input[i])[2] >= minWeight) {
        output[j++] = [point[0], point[1]];
      }
    }

    output.length = j;
    return output;
  });

  return {
    type: "Topology",
    transform: topology.transform,
    bbox: topology.bbox,
    objects: topology.objects,
    arcs: arcs
  };
};

var pi = Math.PI;
var tau = 2 * pi;
var quarterPi = pi / 4;
var radians = pi / 180;
var atan2 = Math.atan2;
var cos = Math.cos;
var sin = Math.sin;

function halfArea(ring, closed) {
  var i = 0,
      n = ring.length,
      sum = 0,
      point = ring[closed ? i++ : n - 1],
      lambda0, lambda1 = point[0] * radians,
      phi1 = (point[1] * radians) / 2 + quarterPi,
      cosPhi0, cosPhi1 = cos(phi1),
      sinPhi0, sinPhi1 = sin(phi1);

  for (; i < n; ++i) {
    point = ring[i];
    lambda0 = lambda1, lambda1 = point[0] * radians;
    phi1 = (point[1] * radians) / 2 + quarterPi;
    cosPhi0 = cosPhi1, cosPhi1 = cos(phi1);
    sinPhi0 = sinPhi1, sinPhi1 = sin(phi1);

    // Spherical excess E for a spherical triangle with vertices: south pole,
    // previous point, current point.  Uses a formula derived from Cagnoli’s
    // theorem.  See Todhunter, Spherical Trig. (1871), Sec. 103, Eq. (2).
    // See https://github.com/d3/d3-geo/blob/master/README.md#geoArea
    var dLambda = lambda1 - lambda0,
        sdLambda = dLambda >= 0 ? 1 : -1,
        adLambda = sdLambda * dLambda,
        k = sinPhi0 * sinPhi1,
        u = cosPhi0 * cosPhi1 + k * cos(adLambda),
        v = k * sdLambda * sin(adLambda);
    sum += atan2(v, u);
  }

  return sum;
}

function sphericalRingArea(ring, interior) {
  var sum = halfArea(ring, true);
  if (interior) sum *= -1;
  return (sum < 0 ? tau + sum : sum) * 2;
}

function sphericalTriangleArea(t) {
  var sum = halfArea(t, false);
  return (sum < 0 ? tau + sum : sum) * 2;
}

exports.bbox = bbox;
exports.feature = feature;
exports.mesh = mesh;
exports.meshArcs = meshArcs;
exports.merge = merge;
exports.mergeArcs = mergeArcs;
exports.neighbors = neighbors;
exports.quantize = quantize;
exports.transform = transform;
exports.untransform = untransform;
exports.topology = topology;
exports.filter = filter;
exports.filterAttached = filterAttached;
exports.filterAttachedWeight = filterAttachedWeight;
exports.filterWeight = filterWeight;
exports.planarRingArea = planarRingArea$1;
exports.planarTriangleArea = planarTriangleArea;
exports.presimplify = presimplify;
exports.quantile = quantile;
exports.simplify = simplify;
exports.sphericalRingArea = sphericalRingArea;
exports.sphericalTriangleArea = sphericalTriangleArea;

Object.defineProperty(exports, ''__esModule'', { value: true });

})));');

create or replace function plv8.plv8_startup()
returns void
language plv8
as
$$
load_module = function(modname) {
 var rows = plv8.execute("SELECT code from plv8_modules " +" where modname = $1", [modname]);
 for (var r = 0; r < rows.length; r++) {
	    var code = rows[r].code;
	    eval("(function() { " + code + "})")();
	  }
	};
	$$;

	select plv8.plv8_startup();

