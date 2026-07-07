function [EQ_real,EQ_imag] = evaluateBeltramiEquationResidual(v_k)
% EVALUATEBELTRAMIEQUATIONRESIDUAL Nonlinear equation used in the Beltrami-based reconstruction step.
% Author: Joey.Chow
    global A_2F A_3F B_2F B_3F Bel_coeff
    UPP =   (A_2F+A_3F*v_k(1)-B_3F*v_k(2))+sqrt(-1)*(A_3F *v_k(2)+B_2F+B_3F*v_k(1));
    DOWNN = (A_2F+A_3F*v_k(1)+B_3F*v_k(2))+sqrt(-1)*(A_3F *v_k(2)-B_2F-B_3F*v_k(1));
    EQ_real = real(UPP./DOWNN)-real(Bel_coeff);
    EQ_imag = imag(UPP./DOWNN)-imag(Bel_coeff);
end
