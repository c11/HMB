AutoCorrSum = function(N, phi) {
  
  AutoCorrSum_seq = rep(0, length(phi));
  for(k in 1:length(phi)) {
    
CorrSum2by2 = 4*phi[k]^1 + 2*phi[k]^sqrt(2);
CorrSum3by3 = 4*CorrSum2by2 + 6*phi[k]^2 + 8*phi[k]^sqrt(5) + 2*phi[k]^sqrt(8) - 4*phi[k]^1;
CorrSum2by1 = phi^1;
CorrSum3by2 = 7*phi[k]^1 + 2*phi[k]^2 + 4*phi[k]^sqrt(2) + 2*phi[k]^sqrt(5);
CorrSumQuadPrev2 = CorrSum2by2;
CorrSumQuadPrev1 = CorrSum3by3;
CorrSumRectPrev2 = CorrSum2by1;
CorrSumRectPrev1 = CorrSum3by2;

for (i in 4:sqrt(N)) {
  CorrSumQuad0 = 0;
  CorrSumQuad0 = 4*CorrSumQuadPrev1 - 4*CorrSumRectPrev1 + CorrSumQuadPrev2;
  ExtraQuad = 0;
  ExtraRect = 0;
  
  for (j in 0:(i-1)) {
    if (j == 0) {
      ExtraQuad = ExtraQuad + 2*i*phi[k]^(i-1);
    } else {
      if (j == (i-1)){
        ExtraQuad = ExtraQuad + 2*phi[k]^sqrt((i-1)^2+(i-1)^2);
      } else {
        ExtraQuad = ExtraQuad + 4*(i-j)*phi[k]^sqrt((i-1)^2 + j^2); 
      } } }
  
  for (j in 0:(i-2)){
    if (j==0) {
      ExtraRect = ExtraRect + (i-1)*phi[k]^(i-1);
    } else {
      ExtraRect = ExtraRect + 2 * (i-1-j)*phi[k]^sqrt((i-1)^2 + j^2);
    } }
  
  CorrSumQuad0 =CorrSumQuad0 + ExtraQuad;
  CorrSumRectPrev2 = CorrSumRectPrev1;
  CorrSumRectPrev1 = 2*CorrSumQuadPrev1 - CorrSumRectPrev2 + ExtraRect; 
  CorrSumQuadPrev2 = CorrSumQuadPrev1;
  CorrSumQuadPrev1 = CorrSumQuad0;
  } 
AutoCorrSum_seq[k] = 2*CorrSumQuad0;   
  }
return(AutoCorrSum_seq); 
}

AutoCorrSum2 = function(N, phi) {
    
  AutoCorrSum_seq2 = rep(0, length(N));
  for(k in 1:length(N)) {
    
    CorrSum2by2 = 4*phi^1 + 2*phi^sqrt(2);
    CorrSum3by3 = 4*CorrSum2by2 + 6*phi^2 + 8*phi^sqrt(5) + 2*phi^sqrt(8) - 4*phi^1;
    CorrSum2by1 = phi^1;
    CorrSum3by2 = 7*phi^1 + 2*phi^2 + 4*phi^sqrt(2) + 2*phi^sqrt(5);
    CorrSumQuadPrev2 = CorrSum2by2;
    CorrSumQuadPrev1 = CorrSum3by3;
    CorrSumRectPrev2 = CorrSum2by1;
    CorrSumRectPrev1 = CorrSum3by2;
    
    for (i in 4:sqrt(N[k])) {
      CorrSumQuad0 = 0;
      CorrSumQuad0 = 4*CorrSumQuadPrev1 - 4*CorrSumRectPrev1 + CorrSumQuadPrev2;
      ExtraQuad = 0;
      ExtraRect = 0;
      
      for (j in 0:(i-1)) {
        if (j == 0) {
          ExtraQuad = ExtraQuad + 2*i*phi^(i-1);
        } else {
          if (j == (i-1)){
            ExtraQuad = ExtraQuad + 2*phi^sqrt((i-1)^2+(i-1)^2);
          } else {
            ExtraQuad = ExtraQuad + 4*(i-j)*phi^sqrt((i-1)^2 + j^2); 
          } } }
      
      for (j in 0:(i-2)){
        if (j==0) {
          ExtraRect = ExtraRect + (i-1)*phi^(i-1);
        } else {
          ExtraRect = ExtraRect + 2 * (i-1-j)*phi^sqrt((i-1)^2 + j^2);
        } }
      
      CorrSumQuad0 =CorrSumQuad0 + ExtraQuad;
      CorrSumRectPrev2 = CorrSumRectPrev1;
      CorrSumRectPrev1 = 2*CorrSumQuadPrev1 - CorrSumRectPrev2 + ExtraRect; 
      CorrSumQuadPrev2 = CorrSumQuadPrev1;
      CorrSumQuadPrev1 = CorrSumQuad0;
    } 
    AutoCorrSum_seq2[k] = 2*CorrSumQuad0;   
  }
  return(AutoCorrSum_seq2); 
}