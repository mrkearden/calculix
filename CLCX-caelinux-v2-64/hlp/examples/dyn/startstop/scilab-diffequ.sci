debug

// diff equation
// y'' = -omega^2*y-2*eta*omega*y'+u(t)/m
//or
//m*y''+C*y'+K*y=u(t)

//Static Force Value
fmax=220

// damping coefficient C/(2*m*omega)
//https://en.wikipedia.org/wiki/Damping_ratio
eta=0.01

//Frequency Of Structure (exact value for given input)
omega=0.73*2*%pi

//Frequency of periodic force (make sense p>omega)
p=1.5*2*%pi

//mass
m=251.2

//start time
ts=3.3
//process time
tp=3.35
//stopping time
te=10

//Full time
Tm=ts+tp+te
//Initial BC
y0=0
v0=0
//Acceleration of start and stop
epss=p/ts
epse=p/te

//Control K value (rigidity of spring or structure)
Keq=omega^2*m

//Amplitude of force at end of start time (for smooth transition)
samp=(epss*sin(epss*ts^2/2)+(epss*ts)^2*cos(epss*ts^2/2))/p^2


//unbalanced rotating force per time
function ft=u(t)
    ft=fmax*cos(p*t-acos(samp));
    if t<=ts then ft=fmax*(epss*sin(epss*t^2/2)+(epss*t)^2*cos(epss*t^2/2))/p^2;end
    if t>=ts+tp then ft=fmax*(epse*sin(epse*(Tm-t)^2/2)+(epse*(Tm-t))^2*cos(epse*(Tm-t)^2/2))/p^2;end      
endfunction

//plot it
//x=[0:ts/1000:Tm];
//plot(x,u)

//Solve ODE
function dz=f(t,y)
    dz=[y(2),-omega^2*y(1)-2*eta*omega*y(2)+u(t)/m]
endfunction


z0=[y0;v0];t0=0;t=0:Tm/5000:Tm;
y=ode(z0,t0,t,f)


y1=y(1,:)
plot(t,abs(y1))


MaxY=max(abs(y1))
StaticY=fmax/Keq
DynRatio=MaxY/StaticY

