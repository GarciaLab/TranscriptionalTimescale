function [ET12, ET21] = pt_solve(q_in,s1,s2)

  dim = size(q_in,1);
  
  % calculate ss matrix
  [V,D] = eig(q_in');
  [~,mi] = max(diag(D));
  ss = V(:,mi)/sum(V(:,mi));
  PI = repmat(ss',dim,1);
  
  % calculate fundamental matrix
  Z = inv(PI-q_in) - PI;
  
  % calculate passage times from Z
  ET21 = (Z(s1,s1) - Z(s2,s1))/ss(s1);
  ET12 = (Z(s2,s2) - Z(s1,s2))/ss(s2);
end