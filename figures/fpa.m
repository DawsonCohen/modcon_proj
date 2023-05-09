mu = rand(1,2);
Sigma = 0.5*[1 0; 0 1];

x1 = linspace(-3,3,64);
x2 = linspace(-3,3,64);
[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];

y = mvnpdf(X,mu,Sigma);
y = reshape(y,length(x2),length(x1));

imagesc(y)
saveas(gca,'laser.png')