% plot lung approximation from PhysiCell

xdim = [0 8000];
ydim = [0 5000];

broch_rad =  235/2;
broch_cent1 = [1500 2500];
broch_cent2 = [2500 3500];
broch_cent3 = [4000 4000];
broch_cent4 = [5500 3000];
broch_cent5 = [6500 2200];

%X = lhsdesign(12,1)*(7500-500)+500;
%Y = lhsdesign(12,1)*(4500-500)+500;
broch_rad2 =  110/2;

X = 1.0e+03*[4.4214,5.9527,2.2982,4.8693,5.1817,1.5516,3.9351,0.9099,3.1854,6.6712,7.3183,1.8453];
Y = 1.0e+03*[3.8104,3.4388,1.5876,0.6011,2.2410,0.9994,2.6915,1.9629,1.2147,4.2764,2.9921,4.0671];

broch_cent6 = [X(1),Y(1)];
broch_cent7 = [X(2),Y(2)];
broch_cent8 = [X(3),Y(3)];
broch_cent9 = [X(4),Y(4)];
broch_cent10 = [X(5),Y(5)];
broch_cent11 = [X(6),Y(6)];
broch_cent12 = [X(7),Y(7)];
broch_cent13 = [X(8),Y(8)];
broch_cent14 = [X(9),Y(9)];
broch_cent15 = [X(10),Y(10)];
broch_cent16 = [X(11),Y(11)];
broch_cent17 = [X(12),Y(12)];
broch_cent18 = [X(10),Y(10)];
broch_cent19 = [X(11),Y(11)];
broch_cent20 = [X(12),Y(12)];

broch1_x = broch_rad*cos([0:0.1:2*pi])+broch_cent1(1);
broch1_y = broch_rad*sin([0:0.1:2*pi])+broch_cent1(2);
broch2_x = broch_rad*cos([0:0.1:2*pi])+broch_cent2(1);
broch2_y = broch_rad*sin([0:0.1:2*pi])+broch_cent2(2);
broch3_x = broch_rad*cos([0:0.1:2*pi])+broch_cent3(1);
broch3_y = broch_rad*sin([0:0.1:2*pi])+broch_cent3(2);
broch4_x = broch_rad*cos([0:0.1:2*pi])+broch_cent4(1);
broch4_y = broch_rad*sin([0:0.1:2*pi])+broch_cent4(2);
broch5_x = broch_rad*cos([0:0.1:2*pi])+broch_cent5(1);
broch5_y = broch_rad*sin([0:0.1:2*pi])+broch_cent5(2);
broch6_x = broch_rad2*cos([0:0.1:2*pi])+broch_cent6(1);
broch6_y = broch_rad2*sin([0:0.1:2*pi])+broch_cent6(2);
broch7_x = broch_rad2*cos([0:0.1:2*pi])+broch_cent7(1);
broch7_y = broch_rad2*sin([0:0.1:2*pi])+broch_cent7(2);
broch8_x = broch_rad2*cos([0:0.1:2*pi])+broch_cent8(1);
broch8_y = broch_rad2*sin([0:0.1:2*pi])+broch_cent8(2);
broch9_x = broch_rad2*cos([0:0.1:2*pi])+broch_cent9(1);
broch9_y = broch_rad2*sin([0:0.1:2*pi])+broch_cent9(2);
broch10_x = broch_rad2*cos([0:0.1:2*pi])+broch_cent10(1);
broch10_y = broch_rad2*sin([0:0.1:2*pi])+broch_cent10(2);
broch11_x = broch_rad2*cos([0:0.1:2*pi])+broch_cent11(1);
broch11_y = broch_rad2*sin([0:0.1:2*pi])+broch_cent11(2);
broch12_x = broch_rad2*cos([0:0.1:2*pi])+broch_cent12(1);
broch12_y = broch_rad2*sin([0:0.1:2*pi])+broch_cent12(2);
broch13_x = broch_rad2*cos([0:0.1:2*pi])+broch_cent13(1);
broch13_y = broch_rad2*sin([0:0.1:2*pi])+broch_cent13(2);
broch14_x = broch_rad2*cos([0:0.1:2*pi])+broch_cent14(1);
broch14_y = broch_rad2*sin([0:0.1:2*pi])+broch_cent14(2);
broch15_x = broch_rad2*cos([0:0.1:2*pi])+broch_cent15(1);
broch15_y = broch_rad2*sin([0:0.1:2*pi])+broch_cent15(2);
broch16_x = broch_rad2*cos([0:0.1:2*pi])+broch_cent16(1);
broch16_y = broch_rad2*sin([0:0.1:2*pi])+broch_cent16(2);
broch17_x = broch_rad2*cos([0:0.1:2*pi])+broch_cent17(1);
broch17_y = broch_rad2*sin([0:0.1:2*pi])+broch_cent17(2);


figure 
fill([xdim(1) xdim(1) xdim(2) xdim(2)],[ydim,flip(ydim)],'b')
hold on 
fill(broch1_x,broch1_y,[0.95 0.95 0.95])
fill(broch2_x,broch2_y,[0.95 0.95 0.95])
fill(broch3_x,broch3_y,[0.95 0.95 0.95])
fill(broch4_x,broch4_y,[0.95 0.95 0.95])
fill(broch5_x,broch5_y,[0.95 0.95 0.95])
fill(broch6_x,broch6_y,[0.95 0.95 0.95])
fill(broch7_x,broch7_y,[0.95 0.95 0.95])
fill(broch8_x,broch8_y,[0.95 0.95 0.95])
fill(broch9_x,broch9_y,[0.95 0.95 0.95])
fill(broch10_x,broch10_y,[0.95 0.95 0.95])
fill(broch11_x,broch11_y,[0.95 0.95 0.95])
fill(broch12_x,broch12_y,[0.95 0.95 0.95])
fill(broch13_x,broch13_y,[0.95 0.95 0.95])
fill(broch14_x,broch14_y,[0.95 0.95 0.95])
fill(broch15_x,broch15_y,[0.95 0.95 0.95])
fill(broch16_x,broch16_y,[0.95 0.95 0.95])
fill(broch17_x,broch17_y,[0.95 0.95 0.95])

set(gca,'FontSize',18)
xlabel('x dim ({\mu}m)')
ylabel('y dim ({\mu}m)')

(pi*broch_rad^2*5+pi*(broch_rad2)^2*15)/(8000*5000)

