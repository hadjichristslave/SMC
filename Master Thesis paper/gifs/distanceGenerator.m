
mean = [1 1 1];

mean = mean';

title("Distributions infered");

covar = eye(3,3)*.1;

for i=1:.2:5
	mean2 = mean + i-1;
	dist = 0.5* [log(det(covar)/det(covar))+trace(pinv(covar)/covar)+	(mean2-mean)'* pinv(covar) *(mean2-mean)];
	temp = mvnrnd(mean(:),covar,100);
	hold off;
	plot3(temp(:,1),temp(:,2),temp(:,3), "o");
	hold on
	temp = mvnrnd(mean2(:),covar,100);
	plot3(temp(:,1),temp(:,2),temp(:,3), "o");
	#input('prompt');
	path  = strcat ("visualizations/image", num2str(i), ".jpg");
	titleText = strcat("KL Divergence:", num2str(dist));
	title (titleText);
	fontsize=20;
	set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
	saveas( 1,path);

end
view(15,15);
