function [NTCP_conv_flash] = NTCP_conv_flash_same_dose(pp_conv, pp_flash, dose, dose_rate, par, n, NTCP0, model)



%INPUTS

%pp_conv: interpolator for SF(D), CONV-RT

%pp_flash: interpolator for SF(D,R), FLASH-RT

%dose: dose distribution (heterogeneous)

%dose_rate: reference dose rate (typically 100 Gy/s) associated to the
    %voxel with the maximum dose

%par: vector of free parameters [D_50, m, M_D, a, M_R, b, A] (article notation) 

%n: vector of volume parameter values
    %NOTE1: the code uses n(1) to normaliza the dose distribution. In our work this is dose for n=1 and therefore n(1)=1  
    %NOTE2: for n<~0.005 the numerical computation of the gEUD may fail and a different numerical method should be used

%NTCP0: reference NTCP value for CONV-RT

%model: flash effect model, =1 (ROD), =2 (phenomenological model)


%OUTPUTS

%NTCP_conv_flash: Nx3 matrices containing the results of the NTCP calculations using the ROD model, 
    %where N is the number of n values (volume parameter). The first column contains the NTCP values 
    %for conventional RT, the second column contains the NTCP values for FLASH, the third column contains
    %the values of D_50 for each value of n:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


D50=par(1);
m=par(2);

M_D=par(3);
p1=par(4);
M_R=par(5);
p2=par(6);
A=par(7);


%%%%%%%

%minimum and maximum dose for the bisection
a=1;
b=60;

[~, NTCP1]= dicotomia(a, a, dose, n(1), D50, m, NTCP0);
[~, NTCP2]= dicotomia(b, b, dose, n(1), D50, m, NTCP0);

%check that a and b are adequate
while NTCP1>NTCP0 || NTCP2<NTCP0
    if NTCP1>NTCP0
        a=a/2
    elseif NTCP2<NTCP0
        b=b*1.2
    end
    [~, NTCP1]= dicotomia(a, a, dose, n(1), D50, m, NTCP0);
    [~, NTCP2]= dicotomia(b, b, dose, n(1), D50, m, NTCP0);
end


%%%%%%%%

%LOOP n values
for i=1:length(n)    

    if i==1

        %dose normalization for n=1 to reach NTCP_conv=NTCP0
        [factor, NTCP]= dicotomia(a, b, dose, n(i), D50, m, NTCP0);
        NTCP_conv_flash(i,1)=NTCP;
        dose=dose*factor;
        
        if model==1
            d_reducida=dose_efectiva(pp_conv, pp_flash, dose, dose_rate);
        else
            d_reducida=dose.*dose_efectiva_flash_pheno(dose, dose_rate*dose/max(dose), M_D, p1, M_R, p2, A);
        end
        gEUD=g_EUD(d_reducida,1/n(i));
        x=(gEUD-D50)/(m*D50);

        NTCP_conv_flash(i,2)=(erf(x/sqrt(2))+1)/2;
        NTCP_conv_flash(i,3)=D50;


    else

        %re-calculation of D50 to still have NTCP_conv=NTCP0 with the same dose distribution

        %check before the bisection algorithm
        aa=D50/2;
        bb=D50*3;
        [D50_new, NTCP1]= dicotomia_D50(aa, aa, dose, n(i), m, NTCP0);
        [D50_new, NTCP2]= dicotomia_D50(bb, bb, dose, n(i), m, NTCP0);

        while NTCP2>NTCP0 || NTCP1<NTCP0
            if NTCP1<NTCP0
                aa=aa/2
            elseif NTCP2>NTCP0
                bb=bb*1.2
            end                    
            [D50_new, NTCP1]= dicotomia_D50(aa, aa, dose, n(i), m, NTCP0);       
            [D50_new, NTCP2]= dicotomia_D50(bb, bb, dose, n(i), m, NTCP0);
        end


        %%
        [D50_new, NTCP]= dicotomia_D50(aa, bb, dose, n(i), m, NTCP0);
        NTCP_conv_flash(i,1)=NTCP;


        if model==1
            d_reducida=dose_efectiva(pp_conv, pp_flash, dose, dose_rate);
        else
            d_reducida=dose.*dose_efectiva_flash_pheno(dose, dose_rate*dose/max(dose), M_D, p1, M_R, p2, A);
        end
        gEUD=g_EUD(d_reducida,1/n(i));
        x=(gEUD-D50_new)/(m*D50_new);

        NTCP_conv_flash(i,2)=(erf(x/sqrt(2))+1)/2;
        NTCP_conv_flash(i,3)=D50_new;


    end

end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gEUD = g_EUD(d,a)

gEUD=((1/length(d))*sum(d.^a))^(1/a);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [factor, NTCP]= dicotomia(a, b, dose, n, D50, m, NTCP0)

%dose escalation

if a==b
    factor=(a+b)/2;
    d_reducida=factor*dose;
    gEUD=g_EUD(d_reducida,1/n);
    x=(gEUD-D50)/(m*D50);
    %keyboard
    NTCP=(erf(x/sqrt(2))+1)/2;
    return
end

%bisection
for i=1:20
    factor=(a+b)/2;
    d_reducida=factor*dose;
    gEUD=g_EUD(d_reducida,1/n);
    x=(gEUD-D50)/(m*D50);
    NTCP=(erf(x/sqrt(2))+1)/2;

    if NTCP==NTCP0
        return
    elseif NTCP<NTCP0
        a=factor;
    elseif NTCP>NTCP0
        b=factor;
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D50, NTCP]= dicotomia_D50(a, b, dose, n, m, NTCP0)

%renomalization of D_50

if a==b
    D50=(a+b)/2;
    d_reducida=dose;
    gEUD=g_EUD(d_reducida,1/n);
    x=(gEUD-D50)/(m*D50);
    NTCP=(erf(x/sqrt(2))+1)/2;
    return
end

%bisection
for i=1:20
    D50=(a+b)/2;
    d_reducida=dose;
    gEUD=g_EUD(d_reducida,1/n);
    x=(gEUD-D50)/(m*D50);
    NTCP=(erf(x/sqrt(2))+1)/2;

    if NTCP==NTCP0
        return
    elseif NTCP<NTCP0
        b=D50;
    elseif NTCP>NTCP0
        a=D50;
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dose_F]=dose_efectiva(pp_conv, pp_flash, dose, dose_rate)

%effective flash dose (ROD)
%D_F / SF_c(D_F)=SF_f(D)
    maximo=max(dose);

    for i=1:length(dose)

        SF=pp_flash(dose(i),dose_rate*dose(i)/maximo);

        a=dose(i)/2;
        SF1=ppval(pp_conv,a);
        ind=0;
        while SF1<SF
            a=a/2;
            SF1=ppval(pp_conv,a);
            ind=ind+1;
            if ind>10, keyboard, end
        end
        b=dose(i);


        for j=1:20
            c=(a+b)/2;
            SFc=ppval(pp_conv,c);
            if SFc<SF
                b=c;
            else
                a=c;
            end
        end

        dose_F(i)=(a+b)/2;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FMF]=dose_efectiva_flash_pheno(dose, dose_rate, k, n, k2, n2, maximo)

%effective flash dose (phenmenological)

F1=dose_rate.^n2./(k2^n2+dose_rate.^n2);
F2=dose.^n./(k^n+dose.^n);

FMF=(1-maximo*F1.*F2);


end