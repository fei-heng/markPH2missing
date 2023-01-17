function [wei]=weight(Z_c,sel,method)

if ( method=='parametric' )
  theta1 = glmfit(Z_c, sel, 'binomial', 'link', 'logit','constant','on');
  pihat = glmval(theta1,Z_c,'logit','constant','on');
  wei=sel./pihat;
end

if ( method=='nonparamet' )
 kerwei=ksrmv(Z_c,sel,repmat(0.2,size(Z_c,2)));
 kerwei1=kerwei.f;
 wei=sel./kerwei1;
 
end
wei(isnan(wei))=0;