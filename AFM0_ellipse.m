file = 'La_0.csv';
minH = 0;
maxH = 6;
thresh = 2.5;
f = csvread(file,1, 0);
%convert unit to nm
f_nm = f * 1e9;
%find 0 of image
%subtract off real 0
f_final = f_nm - min(min(f_nm));

%plot raw data
figure(51)
imagesc(f_final)
colormap(flipud(pink))
colorbar
caxis([minH, maxH])
set(gcf,'position',[0,0,500,435])

%substract off background
BG = sum(f_final,2) / length(f_final(1,:));
x = linspace(0,5,length(f_final(1,:)));
y = linspace(0,5,length(f_final(:,1)));
pfit = polyfit(y,BG,3);
figure(54)
plot(y, BG)
xlabel('y cross-section (um)')
ylabel('averaged line cuts (nm)')
hold on 
plot(y, pfit(4) + pfit(3)*y + pfit(2)*y.^2 + pfit(1)*y.^3)


[X,Y] = meshgrid(x, y);
f_subtracted = f_final - (pfit(4) + pfit(3) * Y + pfit(2) .* Y.^2 + pfit(1) .* Y.^3);
f_new = f_subtracted - min(min(f_subtracted));

figure(59)
imagesc(f_final - f_new)
colormap(flipud(pink))
colorbar
title('subtracted BG')

figure(55)
imagesc(f_new)
colormap(flipud(pink))
colorbar
caxis([minH, maxH])
set(gcf,'position',[0,0,500,435])

f_filter = f_new>thresh;
figure(53)
imagesc(f_filter)
colormap(flipud(pink))
colorbar
set(gcf,'position',[0,0,500,435])

drop_small = bwareaopen(f_filter,3);
figure(66)
imagesc(drop_small)
colormap(flipud(pink))
colorbar
set(gcf,'position',[0,0,500,435])

cc = bwconncomp(drop_small,4);
labeled = labelmatrix(cc);

%check code identified clusters
identified = labeled > 0;
figure(60)
imagesc(identified)
colormap(flipud(autumn))
colorbar
set(gcf,'position',[0,0,500,435])


% get grain area, center and height 
graindata = regionprops(cc,'all');
grain_Centroid = cat(1,graindata.Centroid);
grain_areas = [graindata.Area];
[junk_area, junk_i] = max(grain_areas);
grain_only_Centroid = grain_Centroid;
grain_only_Centroid(junk_i, :) = [];


figure(61)
imagesc(f_filter)
colormap(flipud(pink))
colorbar
set(gcf,'position',[0,0,500,435])
%hold on
%plot(grain_only_Centroid(:,1),grain_only_Centroid(:,2),'r+')

rounded_centers = round(grain_only_Centroid);
center_only_heights = diag(f_new(rounded_centers(:,2), rounded_centers(:,1)));

% get grain major and minor axes
grain_Major = cat(1,graindata.MajorAxisLength);
grain_Minor = cat(1,graindata.MinorAxisLength);
grain_only_Major = grain_Major;
grain_only_Major(junk_i) = [];
grain_only_Minor = grain_Minor;
grain_only_Minor(junk_i) = [];


% AFM tip angle
alpha = 20;
alpha_rad = alpha/180 * pi;
% AFM diameter
D = 16;
% geain width in nm 
%grain_only_Major_raw_nm = grain_only_Major * 9.78 - 2 * center_only_heights * tan(alpha_rad) - D;
%grain_only_Minor_raw_nm = grain_only_Minor * 9.78 - 2 * center_only_heights * tan(alpha_rad) - D;
grain_only_Major_raw_nm = grain_only_Major * 9.78 - D;
grain_only_Minor_raw_nm = grain_only_Minor * 9.78 - D;
% delete particles with unphysical width and height
to_delete = find(grain_only_Major_raw_nm < 0 | grain_only_Minor_raw_nm < 0);
grain_only_Major_nm = grain_only_Major_raw_nm;
grain_only_Minor_nm = grain_only_Minor_raw_nm;
grain_only_Major_nm(to_delete) = [];
grain_only_Minor_nm(to_delete) = [];
center_only_heights(to_delete) = [];
grain_only_Centroid(to_delete, :) = [];

figure(61)
hold on
plot(grain_only_Centroid(:,1),grain_only_Centroid(:,2),'r+')

figure(65)
plot(grain_only_Major_nm, grain_only_Minor_nm, 'o')
xlabel('major axis length (nm)')
ylabel('minor axis length (nm)')
vols = pi/4 * center_only_heights .* grain_only_Major_nm .* grain_only_Minor_nm;
figure(63)
hg = histogram(vols);
hg.BinWidth = 200;
x_bin = hg.BinEdges + hg.BinWidth/2;
x_bin(end) = [];
y_bin = hg.Values;
xlabel('volume (nm^3)')
set(gca, 'FontSize',20)

% convert to Bohr magneton
Bohr_NiO = 2.8;
Bohr_num = vols/.073 * 4 * Bohr_NiO;
figure(64)
histogram(Bohr_num, 30)
xlabel('m (\mu_B)')
ylabel('counts')
set(gca, 'FontSize',24)

% fit for gamma function
[p,ci] = gamfit(vols);
gamma = gampdf(x_bin,p(1),p(2));
area_hist = trapz(x_bin, y_bin);
area_pdf = trapz(x_bin, gamma);
amplitude = area_hist/area_pdf;

figure(63)
hold on
plot(x_bin,gamma*amplitude,'-r', 'LineWidth',2)

figure(65)
imagesc(f_new)
colormap(flipud(pink))
colorbar
caxis([minH, maxH])
set(gcf,'position',[0,0,500,435])
hold on
plot(grain_only_Centroid(:,1),grain_only_Centroid(:,2),'r+')

vol_frac = sum(vols)/(5000 * 5000 * 7);