clear;
x=[-100,-85,-65,-50,-35,-10,0,20,40,55,65,85,105,125,145,160,175,195];
pin=[57.6,187*6.02/450,0.27,15,21];
ex1=exp(-(x-pin(4))/pin(1));
ex1(1:7)=0
para1=1/(sqrt(2*pi)*pin(5));
para2=2*pin(5)*pin(5);
gx1=para1*exp(-x.*x/para2)
LD1=pin(2)*cov(ex1,gx1)+pin(3)


