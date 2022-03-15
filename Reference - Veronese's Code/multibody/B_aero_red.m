function B_aero_red = B_aero_red(in1,in2,U)
%B_AERO_RED
%    B_AERO_RED = B_AERO_RED(IN1,IN2,U)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    21-Jul-2021 20:00:43

eta1 = in1(1,:);
eta2 = in1(2,:);
eta3 = in1(3,:);
eta4 = in1(4,:);
eta_dot1 = in2(1,:);
eta_dot2 = in2(2,:);
eta_dot3 = in2(3,:);
eta_dot4 = in2(4,:);
t2 = U.^2;
t3 = 1.0./U;
t4 = eta_dot3.*8.450823368442517e-2;
t5 = eta4.*4.428078480963271e-3;
t6 = eta_dot4.*1.091152606114523e-1;
t7 = eta2.*3.67175643271207e-2;
t8 = eta_dot2.*3.67175643271207e-2;
t11 = eta4.*3.914107986960835e-1;
t12 = eta_dot1.*4.677277662084094e-1;
t15 = eta1.*2.881143716827053e-1;
t16 = eta_dot1.*2.881143716827053e-1;
t17 = eta4.*1.450583066977313e-1;
t18 = eta4.*1.150450010913495e-1;
t19 = eta_dot4.*1.150450010913495e-1;
t20 = eta4.*7.022939169265402e-2;
t21 = eta2.*1.549846832782947e-1;
t22 = eta_dot1.*3.29666014863616e-1;
t23 = eta3.*3.113108846761223e-2;
t24 = eta_dot4.*6.429965121275674e-2;
t25 = eta1.*2.065224230329747e-1;
t26 = eta_dot1.*2.065224230329747e-1;
t27 = eta3.*2.761201072216585e-1;
t28 = eta_dot3.*2.761201072216585e-1;
t29 = eta_dot1.*2.607735128233221e-2;
t30 = eta_dot4.*1.921304136011411e-1;
t31 = eta4.*4.081634340443635e-1;
t32 = eta_dot4.*4.081634340443635e-1;
t33 = eta4.*2.822980398361037e-1;
t34 = eta3.*6.702446376573398e-1;
t35 = eta2.*3.083081036198933e-1;
t36 = eta_dot2.*3.083081036198933e-1;
t37 = eta3.*5.337714521681294e-2;
t39 = eta4.*1.793446523041063e-1;
t40 = eta_dot1.*2.030650430480467e-1;
t41 = eta4.*1.839390315059036e-1;
t42 = eta_dot4.*2.152876983903853e-1;
t44 = eta_dot3.*5.681134085626201e-2;
t45 = eta4.*7.822976857340336e-1;
t46 = eta_dot2.*5.223806695666804e-2;
t47 = eta_dot2.*3.604288972478452e-2;
t48 = eta_dot3.*2.366734620102493e-2;
t49 = eta1.*1.288573195802777;
t50 = eta1.*1.565754822909382;
t52 = eta_dot3.*2.095325155414764e-2;
t53 = eta1.*6.120622717051624e-1;
t54 = eta2.*1.421431656800319;
t55 = eta3.*1.16188197812786e-1;
t56 = eta2.*1.836311575613157;
t57 = eta1.*9.911805984601413e-1;
t59 = eta1.*1.085508152754731;
t60 = eta_dot3.*6.281105259597304e-2;
t61 = eta_dot2.*6.6029051994506e-2;
t65 = eta1.*7.558421378911147e-1;
t68 = eta1.*1.703208414811562;
t69 = eta2.*2.120331019224355;
t71 = eta_dot4.*7.19791593164598e-2;
t72 = eta3.*6.427508312284931e-2;
t73 = eta_dot3.*6.427508312284931e-2;
t74 = eta4.*3.371887202988724e-1;
t75 = eta1.*1.318629983958135;
t76 = eta_dot1.*9.794724324367192e-2;
t78 = eta_dot2.*1.925504695338625e-2;
t79 = eta2.*1.027466163216267e-1;
t80 = eta2.*9.000816852493335e-1;
t81 = eta3.*8.794242932387424e-2;
t86 = eta3.*1.371414493669336e-1;
t87 = eta_dot2.*5.533086455330048e-3;
t89 = eta2.*1.910275730030793e-1;
t90 = eta2.*2.102826199564655e-1;
t91 = eta4.*4.091678796153322e-1;
t93 = eta2.*2.261804212631459;
t99 = eta1.*1.412597790328782;
t104 = eta1.*1.386520439046449;
t106 = eta2.*2.158157064117955e-1;
t107 = eta3.*8.083725464839059e-1;
t108 = eta3.*6.690051523887144e-1;
t109 = eta3.*3.604948872444108e-1;
t110 = eta3.*2.502956824689748e-1;
t114 = eta_dot2.*6.91349142113272e-1;
t115 = eta_dot4.*1.412828343435586e-1;
t116 = eta_dot4.*3.00641068835267e-1;
t118 = eta_dot4.*3.401189368913965e-3;
t120 = eta_dot2.*1.410463575982163;
t121 = eta_dot4.*2.168319952059519e-1;
t122 = eta_dot1.*1.202650016457405;
t123 = eta_dot4.*6.008797232215683e-1;
t124 = eta_dot3.*5.14812225206052e-1;
t125 = eta_dot1.*1.308227570583101;
t126 = eta_dot1.*4.701225826476485e-1;
t127 = eta_dot2.*1.091795969861698;
t128 = eta_dot2.*1.62861777454215;
t131 = eta_dot1.*1.012834416158307;
t132 = eta_dot3.*1.922511125216897e-1;
t133 = eta_dot1.*7.613218529548407e-1;
t134 = eta_dot2.*1.737282768505373;
t135 = eta_dot3.*5.138601815291619e-1;
t136 = eta_dot3.*2.768946809128829e-1;
t139 = eta_dot3.*6.209077194641022e-1;
t9 = -t5;
t10 = -t6;
t13 = -t7;
t14 = -t8;
t38 = -t11;
t43 = -t17;
t51 = -t21;
t58 = -t23;
t62 = -t27;
t63 = -t28;
t64 = -t30;
t66 = -t31;
t67 = -t32;
t70 = -t34;
t77 = -t42;
t82 = -t45;
t83 = -t46;
t84 = -t47;
t85 = -t48;
t88 = -t61;
t92 = -t71;
t94 = -t72;
t95 = -t73;
t96 = -t74;
t97 = -t78;
t98 = -t79;
t100 = -t81;
t101 = -t87;
t102 = -t89;
t103 = -t90;
t105 = -t91;
t111 = -t106;
t112 = -t108;
t113 = -t110;
t117 = -t116;
t119 = -t118;
t129 = -t123;
t130 = -t124;
t137 = -t132;
t138 = -t135;
t172 = t121+t122+t128+t136;
t140 = t13+t15+t18+t94;
t141 = t14+t16+t19+t95;
t142 = t4+t40+t77+t84;
t143 = t10+t22+t44+t83;
t144 = t26+t36+t63+t67;
t146 = t12+t24+t85+t88;
t147 = t20+t51+t58+t59;
t149 = t37+t43+t49+t102;
t150 = t60+t64+t76+t97;
t151 = t29+t52+t92+t101;
t153 = t39+t65+t98+t100;
t155 = t55+t96+t103+t104;
t160 = t86+t99+t105+t111;
t170 = t114+t126+t129+t130;
t171 = t115+t120+t131+t137;
t174 = t117+t127+t133+t138;
t175 = t119+t125+t134+t139;
t177 = t3.*t172;
t145 = cos(t140);
t148 = cos(t147);
t152 = cos(t149);
t154 = cos(t153);
t156 = cos(t155);
t163 = cos(t160);
t168 = t3.*t144.*7.68096e-1;
t173 = t3.*t170;
t176 = t3.*t171;
t178 = t3.*t174;
t179 = t3.*t175;
t157 = t145.*(1.27e+2./1.25e+2);
t158 = t148.*(1.27e+2./1.25e+2);
t159 = t148.*(1.27e+2./2.5e+2);
t161 = t152.*(1.27e+2./1.25e+2);
t162 = t152.*(1.27e+2./2.5e+2);
t164 = t154.*(1.27e+2./1.25e+2);
t165 = t154.*(1.27e+2./2.5e+2);
t166 = t156.*(1.27e+2./1.25e+2);
t167 = t156.*(1.27e+2./2.5e+2);
t169 = t163.*(1.27e+2./2.5e+2);
t181 = t3.*t141.*t145.*(1.27e+2./2.5e+2);
t180 = t143.*t159;
t182 = t142.*t162;
t183 = -t181;
t184 = t146.*t165;
t185 = t150.*t167;
t186 = t157+t165;
t187 = t158+t162;
t188 = t151.*t169;
t189 = t159+t164;
t190 = t161+t167;
t191 = t166+t169;
t192 = t141.*t186;
t193 = t143.*t187;
t194 = t146.*t189;
t195 = t142.*t190;
t196 = t157+t189;
t197 = t164+t187;
t198 = t158+t190;
t199 = t150.*t191;
t200 = t161+t191;
t201 = t25+t35+t62+t66+t168+t183;
t202 = t141.*t196;
t203 = t146.*t197;
t204 = t143.*t198;
t205 = t157+t197;
t206 = t142.*t200;
t207 = t164+t198;
t208 = t158+t200;
t209 = t184+t192;
t210 = t3.*t209;
t211 = t141.*t205;
t213 = t146.*t207;
t214 = t157+t207;
t215 = t143.*t208;
t216 = t164+t208;
t222 = t180+t194+t202;
t212 = -t210;
t217 = t141.*t214;
t218 = t146.*t216;
t219 = t157+t216;
t223 = t3.*t222;
t226 = t182+t193+t203+t211;
t220 = t141.*t219;
t221 = t53+t70+t80+t82+t173+t212;
t224 = -t223;
t227 = t3.*t226;
t230 = t185+t195+t204+t213+t217;
t225 = t38+t54+t57+t112+t178+t224;
t228 = -t227;
t231 = t3.*t230;
t234 = t188+t199+t206+t215+t218+t220;
t229 = t41+t56+t75+t113+t176+t228;
t232 = -t231;
t235 = t3.*t234;
t233 = t33+t50+t69+t109+t177+t232;
t236 = -t235;
t237 = t9+t68+t93+t107+t179+t236;
et1 = t2.*t145.*t201.*pi.*3.331383606088268e-1+t2.*t148.*t225.*pi.*3.81183330420801e-1+t2.*t154.*t221.*pi.*5.408201622704801e-1;
et2 = t2.*t152.*t229.*pi.*2.347982682810678e-1+t2.*t156.*t233.*pi.*1.132535799924816e-1+t2.*t163.*t237.*pi.*3.015249119465608e-2;
et3 = t2.*t186.*t221.*pi.*6.557841744268244e-1+t2.*t189.*t225.*pi.*1.064606618642677+t2.*t187.*t229.*pi.*7.503608866551202e-1;
et4 = t2.*t196.*t225.*pi.*6.557841744268244e-1+t2.*t190.*t233.*pi.*4.622013155139129e-1+t2.*t197.*t229.*pi.*1.064606618642677;
et5 = t2.*t191.*t237.*pi.*2.229401180954363e-1+t2.*t198.*t233.*pi.*7.503608866551202e-1+t2.*t205.*t229.*pi.*6.557841744268244e-1;
et6 = t2.*t200.*t237.*pi.*4.622013155139129e-1+t2.*t207.*t233.*pi.*1.064606618642677+t2.*t208.*t237.*pi.*7.503608866551202e-1;
et7 = t2.*t214.*t233.*pi.*6.557841744268244e-1+t2.*t216.*t237.*pi.*1.064606618642677+t2.*t219.*t237.*pi.*6.557841744268244e-1;
et8 = t2.*t145.*t201.*pi.*(-4.245546348155455e-2)-t2.*t148.*t225.*pi.*6.040137423787922e-2-t2.*t154.*t221.*pi.*7.634749355103107e-2;
et9 = t2.*t152.*t229.*pi.*(-4.167535664532143e-2)-t2.*t156.*t233.*pi.*2.226405693694915e-2-t2.*t163.*t237.*pi.*6.39774871371404e-3;
et10 = t2.*t186.*t221.*pi.*(-8.357374701093415e-2)-t2.*t189.*t225.*pi.*1.502903416358879e-1-t2.*t187.*t229.*pi.*1.189003429879512e-1;
et11 = t2.*t196.*t225.*pi.*(-8.357374701093415e-2)-t2.*t190.*t233.*pi.*8.203810363252251e-2-t2.*t197.*t229.*pi.*1.502903416358879e-1;
et12 = t2.*t191.*t237.*pi.*(-4.382688373415186e-2)-t2.*t198.*t233.*pi.*1.189003429879512e-1-t2.*t205.*t229.*pi.*8.357374701093415e-2;
et13 = t2.*t200.*t237.*pi.*(-8.203810363252251e-2)-t2.*t207.*t233.*pi.*1.502903416358879e-1-t2.*t208.*t237.*pi.*1.189003429879512e-1;
et14 = t2.*t214.*t233.*pi.*(-8.357374701093415e-2)-t2.*t216.*t237.*pi.*1.502903416358879e-1-t2.*t219.*t237.*pi.*8.357374701093415e-2;
et15 = t2.*t145.*t201.*pi.*(-7.43194297961757e-2)+t2.*t148.*t225.*pi.*6.568931930159756e-2-t2.*t154.*t221.*pi.*2.736587164091223e-2;
et16 = t2.*t152.*t229.*pi.*9.771443980094581e-2+t2.*t156.*t233.*pi.*7.262661340955699e-2+t2.*t163.*t237.*pi.*2.42276420693807e-2;
et17 = t2.*t186.*t221.*pi.*(-1.462980901499522e-1)-t2.*t189.*t225.*pi.*5.386982606478785e-2+t2.*t187.*t229.*pi.*1.293096836645621e-1;
et18 = t2.*t196.*t225.*pi.*(-1.462980901499522e-1)+t2.*t190.*t233.*pi.*1.923512594506807e-1-t2.*t197.*t229.*pi.*5.386982606478785e-2;
et19 = t2.*t191.*t237.*pi.*1.429657744282618e-1+t2.*t198.*t233.*pi.*1.293096836645621e-1-t2.*t205.*t229.*pi.*1.462980901499522e-1;
et20 = t2.*t200.*t237.*pi.*1.923512594506807e-1-t2.*t207.*t233.*pi.*5.386982606478785e-2+t2.*t208.*t237.*pi.*1.293096836645621e-1;
et21 = t2.*t214.*t233.*pi.*(-1.462980901499522e-1)-t2.*t216.*t237.*pi.*5.386982606478785e-2-t2.*t219.*t237.*pi.*1.462980901499522e-1;
et22 = t2.*t145.*t201.*pi.*1.330232255891089e-1-t2.*t148.*t225.*pi.*1.261668372362076e-1+t2.*t154.*t221.*pi.*7.434783717185519e-2;
et23 = t2.*t152.*t229.*pi.*(-2.48930973079e-1)-t2.*t156.*t233.*pi.*2.221548707770418e-1-t2.*t163.*t237.*pi.*8.322743149756723e-2;
et24 = t2.*t186.*t221.*pi.*2.618567432856474e-1+t2.*t189.*t225.*pi.*1.463540101808173e-1-t2.*t187.*t229.*pi.*2.483599158193063e-1;
et25 = t2.*t196.*t225.*pi.*2.618567432856474e-1-t2.*t190.*t233.*pi.*4.900216005492125e-1+t2.*t197.*t229.*pi.*1.463540101808173e-1;
et26 = t2.*t191.*t237.*pi.*(-4.373127377500823e-1)-t2.*t198.*t233.*pi.*2.483599158193063e-1+t2.*t205.*t229.*pi.*2.618567432856474e-1;
et27 = t2.*t200.*t237.*pi.*(-4.900216005492125e-1)+t2.*t207.*t233.*pi.*1.463540101808173e-1-t2.*t208.*t237.*pi.*2.483599158193063e-1;
et28 = t2.*t214.*t233.*pi.*2.618567432856474e-1+t2.*t216.*t237.*pi.*1.463540101808173e-1+t2.*t219.*t237.*pi.*2.618567432856474e-1;
B_aero_red = [et1+et2+et3+et4+et5+et6+et7;et8+et9+et10+et11+et12+et13+et14;et15+et16+et17+et18+et19+et20+et21;et22+et23+et24+et25+et26+et27+et28];
