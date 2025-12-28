data gosner;
infile 'C:\Users\user\Desktop\Research activity\Modeling Gosner_bufo bufo\Analysis\gosner.data.csv'
delimiter = ',' MISSOVER DSD lrecl=13106 firstobs=2;
input 
                  id
                  day
                  group
                  gosner
                  gosner_cat;
run;

proc print data=gosner; run;

data gosner;
set gosner;
timeclss=day;
time2=day*day;
time3=day*day*day;
time=day;
run;

*Exploring the correlation structure;
**********************************;
proc sgscatter data=Gosnerwide;
  title "Scatterplot Matrix for Gosner Data";
  matrix v2 v3 v4 v5 v6; 
run; 


*LMM under MAR;
* MODEL 1. Simplified mean structure, interaction, no random effects;

ods graphics on;
proc mixed data=gosner method=ml covtest;
class  timeclss group id;
model gosner = time group time*group /s outpm=genfit outp=pred; *influence(effect=id iter=5 est);
*random intercept day / subject=id type=un; *solution g gcorr v vcorr;
repeated timeclss / subject=id type=simple;
run;


* MODEL 2. Simplified mean structure, interaction, random intercepts;

ods graphics on;
proc mixed data=gosner method=ml covtest;
class  timeclss group id;
model gosner = time group time*group /s outpm=genfit outp=pred; *influence(effect=id iter=5 est);
random intercept / subject=id type=un; *solution g gcorr v vcorr;
repeated timeclss / subject=id type=simple;
run;


* MODEL 3. Simplified mean structure, interaction, random slopes;

ods graphics on;
proc mixed data=gosner method=ML covtest;
class  timeclss group id;
model gosner = time group time*group /s outpm=genfit outp=pred; *influence(effect=id iter=5 est);
random day / subject=id type=un; *solution g gcorr v vcorr;
repeated timeclss / subject=id type=simple;
run;

* MODEL 4. Simplified mean structure, interaction, random intercepts and slopes;

ods graphics on;
proc mixed data=gosner method=ml covtest;
class  timeclss group id;
model gosner = time group time*group /s outpm=genfit outp=pred; *influence(effect=id iter=5 est);
random intercept day / subject=id type=un; *solution g gcorr v vcorr;
repeated timeclss / subject=id type=simple;
run;

* MODEL 5. Simplified mean structure, interaction, random intercepts and slopes, TOEPH(1) and UN structure, time^2;

ods graphics on;
proc mixed data=gosner method=ML covtest;
class  timeclss group id;
model gosner = group time time*group/s outpm=genfit outp=pred; *influence(effect=id iter=5 est);
random intercept day / subject=id type=un solution g gcorr v vcorr;
repeated timeclss / subject=id type=simple; 
run;

************************************************
*********** MODEL 6. FINAL MODEL ;**************
Apply in Gosner data and deleted Gosner data
************************************************;

ods graphics on;
proc mixed data=gosner_deleted method=ML covtest empirical nobound plots=studentpanel(marginal conditional);
class  timeclss group(ref="1") id;
model gosner = group time*time*group time*time*time*group /s outpm=genfit outp=pred influence(effect=id iter=5 est);
random intercept day / subject=id type=un;* solution g gcorr v vcorr;
repeated timeclss / subject=id type=simple;* ods output solutionr=out; 
*estimate "g1 - g2" time*time*time*group -1 1 0 0;
*estimate "g1 - g3" time*time*time*group -1 0 1 0;
*estimate "g1 - g4" time*time*time*group -1 0 0 1;
*estimate "g3 - g4" time*time*time*group 0 0 1 -1;
*estimate "g2 - g3" time*time*time*group 0 -1 1 0;
*estimate "g2 - g4" time*time*time*group 0 -1 0 1;
run;



*******************************************************************
***************Program for sensitivity analysis;*****************
*Test influence of some observations;
data gosner_deleted;
set gosner;
if id=3 then delete;
if id=13 then delete;
if id=20 then delete;
if id=28 then delete;
if id=31 then delete;
run;

*fit model without influences;
proc mixed data=gosner_deleted method=ML covtest;
class  timeclss group(ref="1") id;
model gosner = group time*time*group time*time*time*group /s outpm=genfit outp=pred;*influence(effect=id iter=5 est);
random intercept day / subject=id type=un; *solution g gcorr v vcorr;
repeated timeclss / subject=id type=simple; 
run;


**********************************************************


**************************************************************************
***** Fit LMM with multiple imputed dataset under MAR;


*Create wide format of data;
proc sort data=gosner_deleted; by id group; run;
PROC TRANSPOSE DATA=gosner_deleted OUT=gosner_wide;
BY id group;
ID time;
VAR gosner;
RUN;


*imputation under MCMC; 
proc mi data=gosner_wide seed=4860495 out=gosner_mi simple nimpute=25;
title ’Model multiple imputation’;
var id _1 _8 _15 _22 _29;
by group;
run; 
proc print data=gosner_mi; run;

*long format;
data gosner_long;
set gosner_mi;
array y (5) _1 _8 _15 _22 _29;
do j=1 to 5;
gosner=y(j);
time=j;
timeclss=time;
time2=time*time;
output;
end;
run;

*rescale time again;
data gosner_long;
set gosner_long;
if time=1 then timec=1;
if time=2 then timec=8;
if time=3 then timec=15;
if time=4 then timec=22;
if time=5 then timec=29;
run;


proc sort data=gosner_long; by _imputation_; run;
proc print data=gosner_long; run;

*Fit mixed model under MAR imputation;
proc mixed data=gosner_long method=ml covtest namelen=30;
by _imputation_;
class timeclss group (ref="1") id;
model gosner = group timec*timec*group timec*timec*timec*group / covb  s outpm=genfit outp=pred;
random intercept timec/ subject=id type=un g gcorr v vcorr;
repeated timeclss / subject=id type=simple;
ods output solutionF=mixparms CovB=mixcovb;
run;

*make mixed models all together MI_MAR;
proc mianalyze parms=mixparms covb(effectvar=rowcol)=mixcovb;
title 'Multiple Imputation Analysis After  LMM for Gosner Data';
class group;
modeleffects intercept group timec*timec*group timec*timec*timec*group;
run; 


