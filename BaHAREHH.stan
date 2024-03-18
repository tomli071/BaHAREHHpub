data {
  int<lower=1> NrOfYears;// Number of years in analysis
  int<lower=1> NrOfLans; //Number of counties,
  int<lower=1> NrOfKretsar; //Number of Hunting Management Precincts (HMPs)
  int<lower=0> NrOfNonFullyCovered; //Number of HMPs where not all area are covered by reports
  int<lower=1> Reports; //Number of reports
  real<lower=0> A[Reports]; // Reported areas per team
  int<lower=1> KretsNrSeqID[Reports]; //Sequential HMP IDs in reports
  int<lower=1> LanSeqID[Reports]; //Sequential County ID in  report
  int<lower=1> KretsListLanSeqID[NrOfKretsar]; //Sequential County ID per circuit
  int<lower=1> KretsYearSeq[NrOfKretsar]; //Sequential Year for HMPs
  real<lower=0>KretsListA[NrOfKretsar]; // Estimated Huntable Land (EJL) per HMP
  int<lower=0, upper=1> FirstKretsAppearance[NrOfKretsar]; // Boolean variable indicating whether the its the first time an HMP occurs in the list of combinations of HMPs and years
  int<lower=1> maxFwdconnections; // used for HMPs that changes over time. Indicates the maximum number of HMPs that a defunct HMP is merged into after the change. 
  int<lower=0> FwdConnectionSeqID[NrOfKretsar,maxFwdconnections]; // Indicates which HMP(s) a changed HMP will belong to after the change. 
  real <lower=0 , upper=1> FwdConnectionprop[NrOfKretsar,maxFwdconnections]; //The proportion of a defunct HMP that will belong to each of the resulting HMPs after the change.
  int<lower=0, upper=1> Hasreports[NrOfKretsar]; // Boolean variable, indcating whether there are any reports for each HMP and year combination
  int<lower=0, upper=1> FullyCovered[NrOfKretsar]; // Boolean variable, indcating whether the area of HMP and year combination is fully covered by the reports
  real sumlogxobsnormalized[NrOfKretsar]; //Sum of log-proportion of EHL covered by reporting teams. Used in the Dirichlet model for nonreporting teams. Becasue it is constant, it is precalculated before running the analysis 
  int<lower=0> NonFullyCoveredIndex[NrOfKretsar]; //Indices used for HMPs not fully covered 
  int<lower=0, upper=1> iszeroK[Reports]; // Boolean variable, indicating whether a report is a null-harvest. This is used to simplify calculations for such reports. 
  int numberzeroK; // Number of reports with null-harvest
  int numberNOTzeroK; // Number of reports with harvest >0
  int whichiszeroK[numberzeroK]; // report indices for null-harvest
  int whichisNOTzeroK[numberNOTzeroK]; // report indices for harvest>0
  real<lower=0> M[NrOfKretsar]; //Number of reporting teams per HMP
  real<lower=0, upper=1> Ahat[NrOfKretsar]; // EHL not covered by reports per HMP
  real<upper=0> logAhat[NrOfKretsar]; //log Ahat per HMP
  real logAtot[NrOfKretsar]; //log total reported are per HMP. 
  vector[Reports] xobsall; //proportion of EHL covered by each team. Becasue it is constant, it is precalculated before running the analysis 
  vector[Reports] xobsnormalizedall;//xobsall normalized to sum to one. Becasue it is constant, it is precalculated before running the analysis 
  int reportedKforKretsarSeqID[NrOfKretsar]; // reported harvest per HMP.
  
  
  int<lower=0> K[Reports]; //Reported harvest
  vector<lower=0>[Reports] Aopt; //Reported area in vector format. Not implemented. 
  
  int<lower=0> NrOfNon0LanKretsar; // Number of HMPs for which harvest is analyzed. 
  int<lower=0> NrOfNon0LanLans;// Number of counties for which harvest is analyzed. 
  int<lower=0> NrOfNon0LanAndNotFullycoveredKretsar; // Number of HMPs for which harvest is analyzed and EHL is not covered by reported area.
  int<lower=0> WhichisNon0LanKrets[NrOfNon0LanKretsar]; //Which HMPs are included for harvest parameters. HMPs (and reports) in counties with only null-harvest across all years are excluded from the analysis of harvest parameters (if no harvest has been reported in any year, the harvest rate is assumed to be zero) but is inclued in the analysis of hunting area per team.
  int<lower=0> NonZeroLanKretsIndex[NrOfKretsar];//Index for which HMPs are included for harvest parameters
  int<lower=0> WhichisNon0LanLan[NrOfNon0LanLans];// which Counties are included for harvest parameters
  int<lower=0> NonZeroLanLansIndex[NrOfLans]; //Index for which Counties are included for harvest parameters
  int<lower=0> WhichisNon0LanAndNotFullycoveredKrets[NrOfNon0LanAndNotFullycoveredKretsar]; //NOT USED
  int<lower=0> NonZeroLanAndNotFullycoveredKretsIndex[NrOfKretsar]; //NOT USED
  int<lower=0> ReportsInNONZeroLan; //Nr Of reports  included for harvest parameters 
  int<lower=0> whichisreportinNONzeroLan[ReportsInNONZeroLan]; //Which reports are included for harvest parameters 
  int<lower=0> NonZeroLanReportsIndex[Reports]; //Index for which reports are included for harvest parameters 
  int<lower=0> ReportszeroKinNON0Lan; //Nr of reports with null-harvest in counties included for harvest parameters 
  int<lower=0> ReportsNOTzeroKinNON0Lan; //Nr of reports with harvest >0 
  int<lower=0> whichiszeroKinNON0Lan[ReportszeroKinNON0Lan]; // NOT USED
  int<lower=0> whichisNOTzeroKinNON0Lan[ReportsNOTzeroKinNON0Lan]; //NOT USED
  int<lower=0, upper=1> IncLansPhi; // Boolean variable, indicating whether county effects are included in the non-linearity of the relationship between area and harvest rate. TRUE in described analysis
  int<lower=0, upper=1> IncYearPhi;// Boolean variable, indicating whether year effects are included in the non-linearity of the relationship between area and harvest rate.  TRUE in described analysis
  int<lower=0, upper=1> AnalyzeHunt; //  Boolean variable, indicating whether harvest analysis should be included. If false, only area per team is analyzed
  int<lower=0, upper=1> AnalyzeHuntTrend; //  Boolean variable, indicating whether a linear trend should be included in the analysis of harvest rates. FALSE in described analyses and and not validated
  int<lower=0, upper=1> AnalyzeHuntTrendLan; //  Boolean variable, indicating whether county effects should be included in the linear trend should be included in the analysis of harvest rates.FALSE in described analyses and and not validated
  int<lower=0, upper=1> Incu; //  Boolean variable, indicating whether year effects should be included for the concentration parameter describing variability in team area
  
  
  int<lower=0, upper=1> doopt; //REDUNDANT
  
  
  
  
  int<lower=0> NrOfFwdConnection[NrOfKretsar]; // The number of HMPs that a defunct HMP is merged into after the change. 
  int<lower=0, upper=1> issimplemerge[NrOfKretsar]; // Boolean variable, indicating whether a the change in HMPs is a merge of two or more HMPs into 1
  
  int<lower=1, upper=NrOfYears> YearSeq[Reports]; // year indices
  real KretsYearSeqmmeanyear[NrOfKretsar]; // Year compared to the average year included in the analysis. Only used for trends in harvest rates.
  
  real lowmean; // mean of prior  for sigma-parameters assumed to model low variability
  real<lower=0> lowsd; // standard deviation of prior for sigma-parameters assumed to model low variability
  real medmean; //  mean of prior  for sigma-parameters assumed to model intermediate variability. Not used.
  real<lower=0>  medsd;// standard deviation of prior for sigma-parameters assumed to model intermediate variability. Not used.
  real highmean; // mean of prior for sigma-parameters assumed to model low variability
  real<lower=0>  highsd; // standard deviation of prior for sigma-parameters assumed to model high variability
  
  
  int<lower=0, upper=1> dopredofexc; //NOT USED
  int<lower=0, upper=1> Incnegrho;// Boolean variable, indicating whether negative temporal autocorrelation is included.
  int<lower=0> ExcludedKperkretsSeqID[NrOfKretsar];//NOT USED
}

parameters {
  real logQ[NrOfNonFullyCovered]; // log of Q, the estimated number of nonreporting teams 
  real logsigmaW[NrOfYears>1 ? 1 :0];//log of sigma_{\omega,m} 
  real logsigmau[Incu && NrOfYears>1? 1 : 0];// log of sigma_{omega,a} 
  real logsigmaL;//Standard deviation of county effects on average log hunting area, sigma_{lambda,m} 
  real logsigmaC; //Standard deviation of HMP effects on average log hunting area, sigma_{chi,m} 
  real W1; //omega_{m,1}
  real W_raw[NrOfYears==1 ? 0 : NrOfYears-1];//omega_{m,t}/sigma_{omega,m}. This parmameterization improves computation by avoiding divergen transitions
  real logsigmav; //log of // sigma_{\lambda,a} 
  // real u; //Nationwide effect on loga
  real L_raw[NrOfLans,NrOfYears];// lambda_{m,l,t}/sigma_{lambda,m}
  real v_raw[NrOfLans,NrOfYears];//lambda_{a,l,t}/sigma_{lambda,a}
  
  real lambda_raw[AnalyzeHunt ? NrOfNon0LanLans : 0, AnalyzeHunt ? NrOfYears : 0];//lambda_{mu,l,t}/sigma_{lambda,mu}
  real xi_raw[AnalyzeHunt  && IncLansPhi ? NrOfNon0LanLans : 0, AnalyzeHunt  && IncLansPhi ? NrOfYears : 0]; //lambda_{phi,l,t}/sigma_{lambda,phi}
  real zeta_raw[AnalyzeHunt  ? NrOfNon0LanLans : 0 , AnalyzeHunt  ? NrOfYears : 0];//lambda_{beta,l,t}/sigma_{lambda,beta}
  real chi_raw[AnalyzeHunt ? NrOfNon0LanKretsar  : 0];// chi_{mu,k,t}/sigma_{chi,mu}
  real omega1[AnalyzeHunt ? 1 : 0]; // omega_{mu,1}
  real eta1[AnalyzeHunt  ? 1 : 0];  // omega_{phi,1}
  real kappa1[AnalyzeHunt  ? 1 : 0]; // omega_{beta,1}
  real u1; // omega_{m,1}
  real logsigmalambda[AnalyzeHunt  ? 1 : 0]; //log sigma_{lambda,mu}
  real logsigmachi[AnalyzeHunt  ? 1 : 0];//log sigma_{chi,mu}
  real logsigmaomega[AnalyzeHunt && NrOfYears>1 ? 1 : 0];//log sigma_{omega,mu}
  real logsigmaeta[AnalyzeHunt && IncYearPhi && NrOfYears>1 ? 1 : 0];//log sigma_{omega,phi}
  real logsigmakappa[AnalyzeHunt && NrOfYears>1 ? 1 : 0];//log sigma_{omega,beta}
  real logsigmaxi[AnalyzeHunt  && IncLansPhi ? 1 : 0]; //log sigma_{lambda,phi}
  real logsigmazeta[AnalyzeHunt  ? 1 : 0]; //log sigma_{lambda,beta}
  
  real omega_raw[AnalyzeHunt  && NrOfYears>1?  NrOfYears-1 : 0]; // parameter on the scale at which omega_{mu,t} is updated for t>1
  real eta_raw[IncYearPhi && AnalyzeHunt && NrOfYears>1 ? NrOfYears-1 : 0]; // parameter on the scale at which omega_{phi,t} is updated for t>1
  real kappa_raw[AnalyzeHunt && NrOfYears>1 ?  NrOfYears-1 : 0];// parameter on the scale at which omega_{beta,t} is updated for t>1
  real u_raw[Incu && NrOfYears>1? NrOfYears -1 : 0]; // parameter on the scale at which omega_{m,t} is updated for t>1
  real logit_rholambda[AnalyzeHunt && NrOfYears>1 ? 1 : 0]; // logit of rho_{lambda,mu} or (rho_{lambda,mu}+1)/2 depending on the inclusion of negative correlation in the analysis. 
  real logit_rhoxi[AnalyzeHunt && IncLansPhi && NrOfYears>1? 1 : 0];// logit of rho_{lambda,phi} or (rho_{lambda,phi}+1)/2 depending on the inclusion of negative correlation in the analysis. 
  real logit_rhozeta[AnalyzeHunt && NrOfYears>1 ? 1 : 0];// logit of rho_{lambda,beta} or (rho_{lambda,beta}+1)/2 depending on the inclusion of negative correlation in the analysis. 
  real logit_rhochi[AnalyzeHunt && NrOfYears>1 ? 1 : 0]; // logit of rho_{chi,phi} or (rho_{chi,phi}+1)/2 depending on the inclusion of negative correlation in the analysis. 
  real psi_raw[AnalyzeHunt && AnalyzeHuntTrendLan  ? NrOfNon0LanLans : 0]; //effects of counties on trends in harvest rates. Not implemented.
  real meanpsi[AnalyzeHunt && AnalyzeHuntTrend ? 1 : 0];//mean of trends in harvest rates. Not implemented.
  real logsigmapsi[AnalyzeHunt && AnalyzeHuntTrendLan ? 1 : 0]; //standard deviation of county effects on trends in harvest rates. Not implemented.
  real logit_rhov[NrOfYears>1 ? 1 : 0]; // logit of rho_{lambda,a} or (rho_{lambda,a}+1)/2 depending on the inclusion of negative correlation in the analysis.  
  real logit_rhoC[NrOfYears>1 ? 1 : 0];// logit of rho_{chi,m} or (rho_{chi,m}+1)/2 depending on the inclusion of negative correlation in the analysis.  
  real logit_rhoL[NrOfYears>1 ? 1 : 0];// logit of rho_{lambda,m} or (rho_{lambda,m}+1)/2 depending on the inclusion of negative correlation in the analysis.  
  // real logbeta;
}

transformed parameters {
  real<lower=0> sigmapsi[AnalyzeHunt  && AnalyzeHuntTrendLan ? 1 : 0]; //Not implemented
  real<lower=0> sigmaomega[AnalyzeHunt && NrOfYears>1 ? 1 : 0];// sigma_{omega,mu}
  real<lower=0> sigmaeta[AnalyzeHunt && IncYearPhi && NrOfYears>1? 1 : 0];//sigma_{omega,phi}
  real<lower=0> sigmakappa[AnalyzeHunt && NrOfYears>1 ? 1 : 0];//sigma_{omega,beta}
  
  
  
  
  real rhoL[NrOfYears>1? 1 : 0];//rho_{lambda,m}
  real rhoC[NrOfYears>1? 1 : 0];//rho_{chi,m}
  real rhov[NrOfYears>1? 1 : 0];//rho_{lambda,a}
  real rholambda[AnalyzeHunt && NrOfYears>1 ? 1 : 0]; //rho_{lambda,mu}
  real rhoxi[AnalyzeHunt && IncLansPhi && NrOfYears>1? 1 : 0]; //rho_{lambda,phi}
  real rhozeta[AnalyzeHunt  && NrOfYears>1? 1 : 0]; //rho_{lambda,beta}
  real rhochi[AnalyzeHunt && NrOfYears>1 ? 1 : 0];////rho_{chi,mu}
  
  
  
  
  real logbetap1[AnalyzeHunt  ? NrOfNon0LanLans : 0 , AnalyzeHunt  ? NrOfYears : 0 ]; // log (beta_{l,t}+1); Calculated once per model evaluation to avoid repeated recalculation for every report
  real logbetamlogbetap1[AnalyzeHunt  ? NrOfNon0LanLans : 0 , AnalyzeHunt  ? NrOfYears : 0 ];// log (beta_{l,t}+1)-log (beta_{l,t}+1); Calculated once per model evaluation to avoid repeated recalculation for every report
  real<lower=0> phi[AnalyzeHunt  ? NrOfNon0LanLans : 0 , AnalyzeHunt  ? NrOfYears : 0 ];//phi_{l,t}
  real<lower=0> beta[AnalyzeHunt  ? NrOfNon0LanLans : 0 , AnalyzeHunt  ? NrOfYears : 0 ];//beta_{l,t}
  real logphi[AnalyzeHunt  ? NrOfNon0LanLans : 0 , AnalyzeHunt  ? NrOfYears : 0 ];// log phi_{l,t}
  real logbeta[AnalyzeHunt  ? NrOfNon0LanLans : 0 , AnalyzeHunt  ? NrOfYears : 0 ];// log beta_{l,t}
  real L[NrOfLans,NrOfYears];//lambda_{m,l}
  real v[NrOfLans,NrOfYears];//lambda_{a,l}
  real C[NrOfKretsar];//chi_{m,l}
  real C_raw[NrOfKretsar];//used in defining C based on Q and unreported area
  real W[ NrOfYears]; //omega_{m,t}
  real u[NrOfYears]; //omega_{a,t}
  real<lower=0> a[NrOfLans,NrOfYears]; //a_{l,t}
  real loga[NrOfLans,NrOfYears]; //log a_{l,t}
  real eta[AnalyzeHunt && IncYearPhi ? NrOfYears : 0];// omega_{phi,t}
  real kappa[AnalyzeHunt  ? NrOfYears : 0]; // omega_{beta,t}
  real xi[AnalyzeHunt &&  IncLansPhi ? NrOfNon0LanLans : 0 , AnalyzeHunt &&  IncLansPhi ? NrOfYears : 0]; //lambda_{phi,l,t}
  real zeta[AnalyzeHunt ? NrOfNon0LanLans : 0 , AnalyzeHunt  ? NrOfYears : 0]; //lambda_{beta,l,t}
  
  vector<lower=0>[NrOfKretsar] Q;//Q
  vector<lower=0>[NrOfKretsar] N;//N
  real logN[NrOfKretsar];// log N
  real<lower=0> sigmaW[NrOfYears>1 ? 1 : 0];//sigma_{omega,m}
  real<lower=0> sigmau[Incu && NrOfYears>1? 1 : 0];//sigma_{omega,a}
  real <lower=0> sigmav;//sigma_{lambda,a}
  real<lower=0> sigmaL;//sigma_{lambda,m}
  real<lower=0> sigmaC; ////sigma_{chi,m}
  
  real sigmalambda[AnalyzeHunt  ? 1 : 0];////sigma_{lambda,mu}
  real<lower=0> sigmachi[AnalyzeHunt  ? 1 : 0];//sigma_{chi,mu}
  real omega[NrOfYears]; //omega_{mu,t}
  real logeffect[AnalyzeHunt  ? NrOfNon0LanKretsar : 0]; // log mu_{k,t} in implemented model
  real logmu[AnalyzeHunt  ? NrOfNon0LanKretsar : 0]; //log mu_{k,t}
  real<lower=0> mu[AnalyzeHunt  ? NrOfNon0LanKretsar : 0];  //mu_{k,t}
  real chi[AnalyzeHunt  ? NrOfNon0LanKretsar : 0];  // chi_{mu,k,t}
  real lambda[AnalyzeHunt ? NrOfNon0LanLans : 0, AnalyzeHunt ? NrOfYears : 0];  // lambda_{mu,l,t}
  
  real nu[AnalyzeHunt  ? ReportsInNONZeroLan : 0]; //nu_i
  
  real alpha[AnalyzeHunt  ? ReportsInNONZeroLan : 0]; //alpha_i
  real alphasum[AnalyzeHunt ? NrOfNon0LanLans : 0, !doopt ? NrOfYears : 0]; // sum of alpha_i for reports i in county and year
  
  real lgammaa [NrOfLans ,  NrOfYears ]; //used for precalculation of the dirichlet parameters, evaluated only once per county and year
  
  real <lower=0> sigmaxi[AnalyzeHunt && IncLansPhi ? 1 : 0]; //sigma_{lambda,phi}
  real <lower=0> sigmazeta[AnalyzeHunt ? 1 : 0 ];//sigma_{lambda,beta}
  real psi[AnalyzeHunt && AnalyzeHuntTrendLan ? NrOfNon0LanLans : 0];// not implemented
  
  
  //////////////////////////////defining rho parameters based on logit rho and sigma parameters based on log sigma
  if (NrOfYears>1){
    if (Incnegrho){
      
      
      rhoL[1]=inv_logit(logit_rhoL[1])*2-1;
      rhoC[1]=inv_logit(logit_rhoC[1])*2-1;
      rhov[1]=inv_logit(logit_rhov[1])*2-1;
      
    }else{
      
      rhoL[1]=inv_logit(logit_rhoL[1]);
      rhoC[1]=inv_logit(logit_rhoC[1]);
      rhov[1]=inv_logit(logit_rhov[1]);
      
    }
    
    
    
    
    
    
  }
  if (AnalyzeHunt){
    if(NrOfYears>1){
      
      if (Incnegrho){
        rhozeta[1]=inv_logit(logit_rhozeta[1])*2-1; 
      }else{
        
        rhozeta[1]=inv_logit(logit_rhozeta[1]); 
        
      }
      
    }
    
    
    
    sigmazeta[1]= exp(logsigmazeta[1]);
    if(AnalyzeHuntTrendLan){
      sigmapsi[1]=exp(logsigmapsi[1]);
    }
    
    if (NrOfYears>1){
      sigmaomega[1]=exp(logsigmaomega[1]);
      
      sigmakappa[1]=exp(logsigmakappa[1]);
      if(Incnegrho){
        
        
        rholambda[1]=inv_logit(logit_rholambda[1])*2-1; 
        
        rhochi[1]=inv_logit(logit_rhochi[1])*2-1;
        
      }else{
        
        
        
        rholambda[1]=inv_logit(logit_rholambda[1]); 
        
        rhochi[1]=inv_logit(logit_rhochi[1]);
        
      }
      
      if (IncYearPhi){
        
        sigmaeta[1]=exp(logsigmaeta[1]);
        
        
      }
    }
    
    if(IncLansPhi){
      if (NrOfYears>1){
        if(Incnegrho){
          rhoxi[1]=inv_logit(logit_rhoxi[1])*2-1; 
        }else{
          rhoxi[1]=inv_logit(logit_rhoxi[1]); 
          
        }
      }
      sigmaxi[1]=exp(logsigmaxi[1]);
    }
    
    
    if(AnalyzeHuntTrendLan){
      
      for (li in 1:NrOfNon0LanLans){
        
        psi[li]=psi_raw[li] * sigmapsi[1];
      }
    }
    
    //////////////////////////////
    
    
    
    ////////////////////////////// defining kappa based on the raw parameters defined for computation
    kappa[1]=kappa1[1];
    if (NrOfYears>1){
      for (iy in 2:NrOfYears){
        
        
        
        kappa[iy]=kappa[iy-1]+kappa_raw[iy-1] * sigmakappa[1];
        
        
      }
      
    }
    
    
    ////////////////////////////// defining xi based on the raw parameters defined for computation
    if (IncLansPhi){
      for (il in 1:NrOfNon0LanLans){
        for (iy in 1:NrOfYears){
          
          
          if (iy == 1){
            xi[il,iy]=xi_raw[il,iy] * sigmaxi[1];
          }else{
            xi[il,iy] = rhoxi[1]*xi[il,iy-1]  + xi_raw[il,iy] * sqrt(1-rhoxi[1]^2) *sigmaxi[1];
            
          }
          
          
          
          
        }
      }
      
      
      
    }
    
    ////////////////////////////// defining eta based on the raw parameters defined for computation
    eta[1]=eta1[1];
    for (iy in 2:NrOfYears){
      
      
      if(IncYearPhi){
        
        eta[iy]=eta[iy-1]+eta_raw[iy-1]*sigmaeta[1];
        
      }else{
        eta[iy]=eta1[1];
        
      }
      
    }
    
    
    ////////////////////////////// defining phi based on underlying parameters
    for (iy in 1:NrOfYears){
      for (il in 1:NrOfNon0LanLans){
        
        
        
        
        
        logphi[il,iy] = eta[iy];
        
        if (IncLansPhi){
          logphi[il,iy] += xi[il,iy];
        }
        
        phi[il,iy] = exp (logphi[il,iy] );
        
      }
    }
    
    ////////////////////////////// defining zeta based on the raw parameters defined for computation
    for (il in WhichisNon0LanLan){
      for (iy in 1:NrOfYears){
        
        
        if (iy == 1){
          zeta[NonZeroLanLansIndex[il],iy]=zeta_raw[NonZeroLanLansIndex[il],iy] * sigmazeta[1];
        }else{
          zeta[NonZeroLanLansIndex[il],iy] = rhozeta[1]*zeta[NonZeroLanLansIndex[il],iy-1] + zeta_raw[NonZeroLanLansIndex[il],iy] * sqrt(1-rhozeta[1]^2) * sigmazeta[1];
        }
        
        
        
        
        
        
      }
    }
    
    
    ////////////////////////////// defining beta based on underlying parameters
    
    for (il in 1:NrOfNon0LanLans){
      for (iy in 1:NrOfYears){
        
        
        
        logbeta[il,iy] =  kappa[iy]+zeta[il,iy];
        
        beta[il,iy] = exp (logbeta[il,iy] );
        
      }
    }
    
    
    
  }
  
  // additional sigma  parameters
  sigmav=exp(logsigmav);
  sigmaC=exp(logsigmaC);
  sigmaL=exp(logsigmaL);
  
  ////////////////////////////// defining W based on the raw parameters defined for computation
  W[1]=W1;
  if (NrOfYears>1){
    sigmaW[1]=exp(logsigmaW[1]);
    for (iy in 2:NrOfYears){
      
      
      
      W[iy]= W[iy-1] + W_raw[iy-1] * sigmaW[1];
      
      
      
      
      
    }
  }
  
  
  ////////////////////////////// defining u based on the raw parameters defined for computation
  u[1]=u1;
  if (NrOfYears>1){
    
    if(Incu){
      sigmau[1]=exp(logsigmau[1]);
      
      for (iy in 2:NrOfYears){
        
        
        
        u[iy]=u[iy-1]+u_raw[iy-1] * sigmau[1];
        
        
        
        
        
        
        
      }
    }else{
      for (iy in 2:NrOfYears){
        u[iy]=u1;
      }
      
    }
  }
  
  
  ////////////////////////////// defining L based on the raw parameters defined for computation
  for (il in 1:NrOfLans){
    for (iy in 1:NrOfYears){
      
      
      if (iy == 1){
        L[il,iy]=L_raw[il,iy] * sigmaL;
      }else{
        L[il,iy] = rhoL[1]*L[il,iy-1]  + L_raw[il,iy] * sqrt(1-rhoL[1]^2) *sigmaL;
        
      }
      
      
      
      
    }
  }
  
  ////////////////////////////// defining v based on the raw parameters defined for computation
  for (il in 1:NrOfLans){
    for (iy in 1:NrOfYears){
      
      
      if (iy == 1){
        v[il,iy]=v_raw[il,iy] * sigmav;
      }else{
        v[il,iy] = rhov[1]*v[il,iy-1]  + v_raw[il,iy] * sqrt(1-rhov[1]^2) *sigmav;
        
      }
      
      
      
      
    }
  }
  
  ////////////////////////////// defining a based on the raw parameters defined for computation
  for (il in 1:NrOfLans){
    for (iy in 1:NrOfYears){
      
      loga[il,iy] = u[iy] +  v[il,iy];
      
      a[il,iy] = exp (loga[il,iy] );
      
    }
  }
  
  
  
  
  ////////////////////////////// defining Q, N and C and based on the raw parameters defined for computation
  for (ik in 1:NrOfKretsar){
    if(FullyCovered[ik]){
      Q[ik]=0;
    }else{
      Q[ik]=exp(logQ[NonFullyCoveredIndex[ik]]);
    }
    
    N[ik]=M[ik]+Q[ik]; 
    logN[ik]=log(N[ik]);
    C[ik]=logAtot[ik]-(W[KretsYearSeq[ik]] +L[KretsListLanSeqID[ik],KretsYearSeq[ik]]+logN[ik]);
    
    if(issimplemerge[ik]){
      real sumA;
      real ACsum;
      sumA=0;
      for (ci in 1:NrOfFwdConnection[ik]){
        sumA += KretsListA[FwdConnectionSeqID[ik,ci]]*FwdConnectionprop[ik,ci];
      }
      
      
      ACsum=0;
      for (ci in 1:NrOfFwdConnection[ik]){
        ACsum += C[FwdConnectionSeqID[ik,ci]]*KretsListA[FwdConnectionSeqID[ik,ci]]*FwdConnectionprop[ik,ci];
        
      }
      
      
      
      C_raw[ik] =(C[ik] - rhoC[1]*ACsum/sumA )/ (sqrt(1-rhoC[1]^2) * sigmaC );
      
      
      
      
    }else{
      if (FirstKretsAppearance[ik]){
        
        C_raw[ik] = C[ik] / sigmaC;
      }else{
        C_raw[ik] =(C[ik] - rhoC[1]*C[FwdConnectionSeqID[ik,1]])/ (sqrt(1-rhoC[1]^2) * sigmaC );
        
      }
      
    }
    
  }
  ////////////////////////////// defining omega based on the raw parameters defined for computation
  if (AnalyzeHunt ){
    omega[1]=omega1[1];
    if(NrOfYears>1){
      
      for (iy in 2:NrOfYears){
        
        
        
        omega[iy] = omega[iy-1] +omega_raw[iy-1] *sigmaomega[1];
        
        
      }
    }
    //additional sigma parameters
    sigmalambda[1]=exp(logsigmalambda[1]);
    sigmachi[1]=exp(logsigmachi[1]);
    
    ////////////////////////////// defining chi based on the raw parameters defined for computation
    for (ik in WhichisNon0LanKrets){
      if(issimplemerge[ik]){
        real sumA;
        real Achisum;
        sumA=0;
        for (ci in 1:NrOfFwdConnection[ik]){
          sumA += KretsListA[FwdConnectionSeqID[ik,ci]]*FwdConnectionprop[ik,ci];
        }
        
        
        Achisum=0;
        for (ci in 1:NrOfFwdConnection[ik]){
          Achisum += chi[NonZeroLanKretsIndex[FwdConnectionSeqID[ik,ci]]]*KretsListA[FwdConnectionSeqID[ik,ci]]*FwdConnectionprop[ik,ci];
          
        }
        
        
        chi[NonZeroLanKretsIndex[ik]] = rhochi[1]*Achisum/sumA + chi_raw[NonZeroLanKretsIndex[ik]] * sqrt(1-rhochi[1]^2) * sigmachi[1]  ;
        
        
        
        
        
        
      }else{
        
        
        if (FirstKretsAppearance[ik]){
          chi[NonZeroLanKretsIndex[ik]]=chi_raw[NonZeroLanKretsIndex[ik]] * sigmachi[1];
        }else{
          chi[NonZeroLanKretsIndex[ik]] = rhochi[1]*chi[NonZeroLanKretsIndex[FwdConnectionSeqID[ik,1]]] + chi_raw[NonZeroLanKretsIndex[ik]] * sqrt(1-rhochi[1]^2) * sigmachi[1] ;
          
        }
        
      }
      
    }
    
    
    ////////////////////////////// defining lambda based on the raw parameters defined for computation
    for (il in WhichisNon0LanLan){
      for (iy in 1:NrOfYears){
        
        
        if (iy == 1){
          lambda[NonZeroLanLansIndex[il],iy]=lambda_raw[NonZeroLanLansIndex[il],iy] * sigmalambda[1];
        }else{
          lambda[NonZeroLanLansIndex[il],iy] = rholambda[1]*lambda[NonZeroLanLansIndex[il],iy-1] + lambda_raw[NonZeroLanLansIndex[il],iy] * sqrt(1-rholambda[1]^2) * sigmalambda[1];
        }
        
        
        
        
        
      }
    }
    
    
    ////////////////////////////// defining mu and logmu based on the raw parameters defined for computation
    for (ik in WhichisNon0LanKrets){
      
      logeffect[NonZeroLanKretsIndex[ik]]=omega[KretsYearSeq[ik]] + lambda[NonZeroLanLansIndex[KretsListLanSeqID[ik]],KretsYearSeq[ik]] + chi[NonZeroLanKretsIndex[ik]] ;
      // if(NrOfYears>1){
        //   logeffect[NonZeroLanKretsIndex[ik]] += omega[KretsYearSeq[ik]] ;
        //   
        // }
        
        if (AnalyzeHuntTrend){
          if (AnalyzeHuntTrendLan){
            
            logeffect[NonZeroLanKretsIndex[ik]] += KretsYearSeqmmeanyear[ik]*(meanpsi[1] + psi[NonZeroLanLansIndex[KretsListLanSeqID[ik]]]);
          }else{
            
            logeffect[NonZeroLanKretsIndex[ik]] += KretsYearSeqmmeanyear[ik]*meanpsi[1];
            
          }
          
          
          
        }
        
        
        
        if (doopt){
          
        }else{
          mu[NonZeroLanKretsIndex[ik]] = exp(logeffect[NonZeroLanKretsIndex[ik]]);
          logmu[NonZeroLanKretsIndex[ik]]=logeffect[NonZeroLanKretsIndex[ik]];
        }
        
        
        
        
        
        
    }
    ////////////////////////////// calculating county level beta parametes once instead of for every report
    for (il in 1:NrOfNon0LanLans){
      for (iy in 1:NrOfYears){
        logbetap1[il,iy]=log(beta[il,iy]+1);
        logbetamlogbetap1[il,iy]=logbeta[il,iy]-logbetap1[il,iy];
      }
    }
    
    ////////////////////////////// defining nu and alpha
    if (doopt){
      
      
    }else{
      for (il in 1:NrOfNon0LanLans){
        for (iy in 1:NrOfYears){
          alphasum[il,iy]=0;
        }
      }
      for (i in whichisreportinNONzeroLan){
        
        nu[NonZeroLanReportsIndex[i]]=KretsListA[KretsNrSeqID[i]]/N[KretsNrSeqID[i]] * mu[NonZeroLanKretsIndex[KretsNrSeqID[i]]]*(N[KretsNrSeqID[i]]*xobsall[i])^phi[NonZeroLanLansIndex[LanSeqID[i]],YearSeq[i]];
        alpha[NonZeroLanReportsIndex[i]]=nu[NonZeroLanReportsIndex[i]]*beta[NonZeroLanLansIndex[LanSeqID[i]],YearSeq[i]];
        if(iszeroK[i]){//This is in order not evaluate the entire negative binomial for reports with K=0
        alphasum[NonZeroLanLansIndex[LanSeqID[i]],YearSeq[i]] += alpha[NonZeroLanReportsIndex[i]];
        }
      }
    }
    
  }
  
  ////////////////////////////// calculating log-gamma of dirichlet concentration parameters once per county
  for (il in 1:NrOfLans){
    for (iy in 1:NrOfYears){
      lgammaa[il,iy]=lgamma(a[il,iy]);
      
    }
  }
  
}
model{
  for (ki in 1:NrOfKretsar){
    
    if (!FullyCovered[ki]){
      
    }
    if(Hasreports[ki]){//Dir process only relevant when there are reports, i.e. M>0
    target += sumlogxobsnormalized[ki]*(a[KretsListLanSeqID[ki],KretsYearSeq[ki]]-1)-M[ki]*lgammaa[KretsListLanSeqID[ki],KretsYearSeq[ki]]+lgamma(M[ki]*a[KretsListLanSeqID[ki],KretsYearSeq[ki]]); //Simplification of Dirichlet likelihood with common shape parameter.
    
    
    if (!FullyCovered[ki]){//Neither Q nor M can be 0 for evaluating the beta  of unreported and total reported areas
    Ahat[ki]~beta(exp(logQ[NonFullyCoveredIndex[ki]] +loga[KretsListLanSeqID[ki],KretsYearSeq[ki]]),M[ki]*a[KretsListLanSeqID[ki],KretsYearSeq[ki]]);
    }
    
    
    
    }
    
    
    
  }
  
  
  // model for Q and necessary hastings correction
  for(ik in 1:NrOfKretsar){
    
    if (FirstKretsAppearance[ik]){
      
      C_raw[ik] ~ normal(0,1);//   reparametrisation for improved computation
      
      if(!FullyCovered[ik]){
        target += logQ[NonFullyCoveredIndex[ik]] - logN[ik] - logsigmaC;
      }
    }else{
      
      C_raw[ik] ~ normal(0, 1);
      if(!FullyCovered[ik]){
        target += logQ[NonFullyCoveredIndex[ik]] - logN[ik] - logsigmaC - log(1-rhoC[1]^2)/2;
        
      }
      //
      
    }
    
    
    
    
    
  }
  
  
  
  /////////////////////////////////////////////////////  defining priors and Hastings corrections for remaining parameters
  for (il in 1:NrOfLans){
    
    
    
    
    
    L_raw[il,1] ~ normal(0,1);
    for (iy in 2:NrOfYears){
      
      
      
      
      L_raw[il,iy] ~ normal(0, 1);
      
      
      
      //
      
      
    }
    
  }
  
  for (il in 1:NrOfLans){
    
    
    
    
    
    v_raw[il] ~ normal(0,1);
    
    
  }
  
  W1 ~ normal(5.8,1.7);
  
  if(NrOfYears>1){
    for (iy in 1:(NrOfYears-1)){
      
      
      
      
      W_raw[iy] ~ normal(0, 1);
      
      
      
      
    }
    
    if(Incu){
      for (iy in 1:(NrOfYears-1)){
        
        
        
        u_raw[iy] ~ normal(0, 1);
        
        
        
        
      }
      
      logsigmau[1] ~ normal(lowmean,lowsd);//
      
    }
  }
  u1 ~ normal(0,1);
  
  
  
  sigmaL ~ cauchy(0,2.5);
  sigmaC ~ cauchy(0,2.5);
  sigmav ~ cauchy(0,2.5);
  
  
  target += logsigmaL;
  target += logsigmaC;
  target += logsigmav;
  
  if(NrOfYears>1){
    logsigmaW[1]~normal(lowmean,lowsd);
    
    if(Incnegrho){
      
      (rhoL[1]+1)/2~beta(1,1);
      (rhoC[1]+1)/2~beta(1,1);
      (rhov[1]+1)/2~beta(1,1);
      
    }else{
      rhoL[1]~beta(1,1);
      rhoC[1]~beta(1,1);
      rhov[1]~beta(1,1);
      
    }
    
    
    target += log_inv_logit(logit_rhoL[1]) + log1m_inv_logit(logit_rhoL[1]);
    target += log_inv_logit(logit_rhoC[1]) + log1m_inv_logit(logit_rhoC[1]);
    target += log_inv_logit(logit_rhov[1]) + log1m_inv_logit(logit_rhov[1]);
    
    
    
    
  }
  
  
  
  
  
  if (AnalyzeHunt){
    if(NrOfYears>1){
      for (iy in 1:(NrOfYears-1)){
        
        
        
        
        
        omega_raw[iy] ~ normal(0, 1);
        if(IncYearPhi){
          eta_raw[iy] ~ normal(0, 1);
        }
        
        kappa_raw[iy] ~ normal(0, 1);
        //
        
        
      }
      
      
      
      
      if(Incnegrho){
        (rhochi[1]+1)/2~beta(1,1);
        (rholambda[1]+1)/2~beta(1,1);
        (rhoxi[1]+1)/2~beta(1,1);
        (rhozeta[1]+1)/2~beta(1,1);
        
        
        
      }else{
        
        rhochi[1]~beta(1,1);
        rholambda[1]~beta(1,1);
        rhoxi[1]~beta(1,1);
        rhozeta[1]~beta(1,1);
        
        
        
      }
      
      
      
      
      if (IncLansPhi){
        target += log_inv_logit(logit_rhoxi[1]) + log1m_inv_logit(logit_rhoxi[1]);
        
      }
      
      
      target += log_inv_logit(logit_rholambda[1]) + log1m_inv_logit(logit_rholambda[1]);
      
      target += log_inv_logit(logit_rhozeta[1]) + log1m_inv_logit(logit_rhozeta[1]);
      target += log_inv_logit(logit_rhochi[1]) + log1m_inv_logit(logit_rhochi[1]);
    }
    if(doopt){
      
      
    }else{
      
      for(i in whichisreportinNONzeroLan){
        if(!iszeroK[i]){
          
          
          K[i]~neg_binomial_2(nu[NonZeroLanReportsIndex[i]],alpha[NonZeroLanReportsIndex[i]]);
        }
      }
      
      for (il in 1:NrOfNon0LanLans){
        for (iy in 1:NrOfYears){
          
          target+= alphasum[il,iy]*logbetamlogbetap1[il,iy];
        }
      }
      
      
    }
    
    
    
    for(ik in WhichisNon0LanKrets){
      
      if (FirstKretsAppearance[ik]){
        chi_raw[NonZeroLanKretsIndex[ik]] ~ normal(0,1);//   
        
        
      }else{
        
        chi_raw[NonZeroLanKretsIndex[ik]] ~ normal(0, 1);
        
        
        
      }
      
      
      
      
    }
    if (AnalyzeHuntTrendLan){
      for (li in 1:NrOfNon0LanLans){
        
        psi_raw[li]~normal(0,1);
      }
      
    }
    
    for (il in WhichisNon0LanLan){
      
      
      
      
      
      lambda_raw[NonZeroLanLansIndex[il] ,1] ~ normal(0,1);
      for (iy in 2:NrOfYears){
        
        lambda_raw[NonZeroLanLansIndex[il] ,iy] ~ normal(0, 1);
        
      }
      
    }
    if (IncLansPhi){
      for (il in WhichisNon0LanLan){
        
        xi_raw[NonZeroLanLansIndex[il]] ~ normal(0,1);
        
        
      }
    }
    for (il in WhichisNon0LanLan){
      
      
      
      
      
      zeta_raw[NonZeroLanLansIndex[il] ] ~ normal(0,1);
      
      
    }
    
    
    sigmalambda[1]~ cauchy(0,2.5);
    target += logsigmalambda[1];
    
    
    sigmachi[1]~ cauchy(0,2.5);
    target += logsigmachi[1];
    
    omega1[1] ~ normal(-8.7,4.3);
    
    eta1[1] ~ normal(0,1);
    kappa1[1] ~ normal(0,5);
    
    
    if(NrOfYears>1){
      logsigmaomega[1]~normal(lowmean,lowsd);
      
    }
    
    if(IncYearPhi && NrOfYears>1){
      
      logsigmaeta[1]~normal(lowmean,lowsd);
      
      
    }
    if(AnalyzeHuntTrendLan){
      
      sigmapsi[1]~ cauchy(0,2.5);
      target += logsigmapsi[1];
      
    }
    
    
    
    if(NrOfYears>1){
      logsigmakappa[1]~normal(lowmean,lowsd);
      
    }
    
    if (IncLansPhi){
      
      sigmaxi[1]~ cauchy(0,2.5);
      target += logsigmaxi[1];
    }
    
    
    sigmazeta[1]~ cauchy(0,2.5);
    target += logsigmazeta[1];
    
  }
}

