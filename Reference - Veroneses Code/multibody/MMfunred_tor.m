function MMred = MMfunred_tor(in1)
%MMFUNRED_TOR
%    MMRED = MMFUNRED_TOR(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    21-Jul-2021 20:00:48

eta1 = in1(1,:);
eta2 = in1(2,:);
eta3 = in1(3,:);
eta4 = in1(4,:);
t2 = eta4.*4.428078480963271e-3;
t3 = eta2.*3.67175643271207e-2;
t5 = eta4.*3.914107986960835e-1;
t7 = eta1.*2.881143716827053e-1;
t8 = eta4.*1.450583066977313e-1;
t9 = eta4.*1.150450010913495e-1;
t10 = eta4.*7.022939169265402e-2;
t11 = eta2.*1.549846832782947e-1;
t12 = eta3.*3.113108846761223e-2;
t13 = eta1.*2.065224230329747e-1;
t14 = eta3.*2.761201072216585e-1;
t15 = eta4.*4.081634340443635e-1;
t16 = eta4.*2.822980398361037e-1;
t17 = eta3.*6.702446376573398e-1;
t18 = eta2.*3.083081036198933e-1;
t19 = eta3.*5.337714521681294e-2;
t21 = eta4.*1.793446523041063e-1;
t22 = eta4.*1.839390315059036e-1;
t24 = eta4.*7.822976857340336e-1;
t25 = eta1.*1.288573195802777;
t26 = eta1.*1.565754822909382;
t28 = eta1.*6.120622717051624e-1;
t29 = eta2.*1.421431656800319;
t30 = eta3.*1.16188197812786e-1;
t31 = eta2.*1.836311575613157;
t32 = eta1.*9.911805984601413e-1;
t34 = eta1.*1.085508152754731;
t36 = eta1.*7.558421378911147e-1;
t38 = eta1.*1.703208414811562;
t39 = eta2.*2.120331019224355;
t41 = eta3.*6.427508312284931e-2;
t42 = eta4.*3.371887202988724e-1;
t43 = eta1.*1.318629983958135;
t44 = eta2.*1.027466163216267e-1;
t45 = eta2.*9.000816852493335e-1;
t46 = eta3.*8.794242932387424e-2;
t48 = eta3.*1.371414493669336e-1;
t49 = eta2.*1.910275730030793e-1;
t50 = eta2.*2.102826199564655e-1;
t51 = eta4.*4.091678796153322e-1;
t52 = eta2.*2.261804212631459;
t56 = eta1.*1.412597790328782;
t60 = eta1.*1.386520439046449;
t62 = eta2.*2.158157064117955e-1;
t63 = eta3.*8.083725464839059e-1;
t64 = eta3.*6.690051523887144e-1;
t65 = eta3.*3.604948872444108e-1;
t66 = eta3.*2.502956824689748e-1;
t4 = -t2;
t6 = -t3;
t20 = -t5;
t23 = -t8;
t27 = -t11;
t33 = -t12;
t35 = -t14;
t37 = -t15;
t40 = -t17;
t47 = -t24;
t53 = -t41;
t54 = -t42;
t55 = -t44;
t57 = -t46;
t58 = -t49;
t59 = -t50;
t61 = -t51;
t67 = -t62;
t68 = -t64;
t69 = -t66;
t83 = t16+t26+t39+t65;
t70 = t6+t7+t9+t53;
t71 = t13+t18+t35+t37;
t74 = t10+t27+t33+t34;
t80 = t28+t40+t45+t47;
t81 = t19+t23+t25+t58;
t82 = t4+t38+t52+t63;
t86 = t21+t36+t55+t57;
t92 = cos(t83);
t93 = sin(t83);
t94 = t22+t31+t43+t69;
t96 = t20+t29+t32+t68;
t103 = t30+t54+t59+t60;
t109 = t48+t56+t61+t67;
t72 = cos(t70);
t73 = sin(t70);
t75 = cos(t71);
t76 = sin(t71);
t77 = cos(t74);
t78 = sin(t74);
t84 = cos(t81);
t85 = sin(t81);
t87 = cos(t80);
t88 = sin(t80);
t89 = sin(t86);
t90 = cos(t82);
t91 = sin(t82);
t95 = cos(t86);
t98 = cos(t94);
t99 = sin(t94);
t100 = cos(t96);
t101 = sin(t96);
t102 = t93.^2;
t105 = cos(t103);
t106 = sin(t103);
t110 = sin(t109);
t111 = cos(t109);
t79 = t76.^2;
t97 = t88.^2;
t104 = t91.^2;
t107 = t99.^2;
t108 = t101.^2;
t112 = t77.*1.574644382107474e-1;
t113 = t78.*1.574644382107474e-1;
t114 = t72.*3.730504535635463e-2;
t115 = t73.*3.730504535635463e-2;
t116 = t72.*2.927242016296286e-1;
t117 = t73.*2.927242016296286e-1;
t118 = t77.*7.873221910537372e-2;
t119 = t78.*7.873221910537372e-2;
t120 = t72.*1.865252267817732e-2;
t121 = t73.*1.865252267817732e-2;
t122 = t72.*1.463621008148143e-1;
t123 = t73.*1.463621008148143e-1;
t124 = t72.*5.844286055440557e-2;
t125 = t72.*1.168857211088111e-1;
t127 = t73.*5.844286055440557e-2;
t128 = t73.*1.168857211088111e-1;
t129 = t72.*6.53034844528149e-2;
t130 = t73.*6.53034844528149e-2;
t135 = t72.*3.265174222640745e-2;
t137 = t73.*3.265174222640745e-2;
t138 = t77.*3.567653097986824e-2;
t139 = t77.*7.135306195973648e-2;
t140 = t78.*3.567653097986824e-2;
t141 = t78.*7.135306195973648e-2;
t142 = t77.*1.581459294154701e-2;
t143 = t77.*3.162918588309402e-2;
t144 = t78.*1.581459294154701e-2;
t145 = t78.*3.162918588309402e-2;
t150 = t77.*5.514381415994032e-1;
t151 = t77.*1.102876283198806;
t152 = t78.*5.514381415994032e-1;
t153 = t78.*1.102876283198806;
t154 = t84.*2.711558977014098e-2;
t155 = t84.*5.423117954028195e-2;
t156 = t85.*2.711558977014098e-2;
t157 = t85.*5.423117954028195e-2;
t161 = t84.*1.940840141711285e-1;
t162 = t84.*7.368961980244748e-2;
t163 = t84.*1.47379239604895e-1;
t164 = t85.*1.940840141711285e-1;
t165 = t85.*7.368961980244748e-2;
t166 = t85.*1.47379239604895e-1;
t167 = t95.*9.110708337048599e-2;
t168 = t95.*1.82214166740972e-1;
t169 = t89.*9.110708337048599e-2;
t170 = t89.*1.82214166740972e-1;
t171 = t84.*6.545951834678109e-1;
t172 = t84.*1.309190366935622;
t173 = t85.*6.545951834678109e-1;
t174 = t85.*1.309190366935622;
t175 = t84.*9.704200708556426e-2;
t178 = t85.*9.704200708556426e-2;
t183 = t95.*3.839678060486863e-1;
t184 = t95.*7.679356120973725e-1;
t185 = t89.*3.839678060486863e-1;
t186 = t89.*7.679356120973725e-1;
t187 = t95.*1.043905621827727e-1;
t188 = t89.*1.043905621827727e-1;
t189 = t89.*4.467475409652812e-2;
t190 = t89.*8.934950819305623e-2;
t191 = t95.*5.219528109138636e-2;
t192 = t89.*5.219528109138636e-2;
t193 = t95.*4.467475409652812e-2;
t194 = t95.*8.934950819305623e-2;
t197 = t105.*1.712918699118272e-1;
t198 = t105.*3.425837398236544e-1;
t199 = t106.*1.712918699118272e-1;
t200 = t106.*3.425837398236544e-1;
t201 = t105.*1.180472089777906e-1;
t202 = t106.*1.180472089777906e-1;
t204 = t105.*1.408704766071192;
t206 = t106.*1.408704766071192;
t207 = t105.*1.068235709378845e-1;
t208 = t105.*2.136471418757689e-1;
t209 = t106.*1.068235709378845e-1;
t210 = t106.*2.136471418757689e-1;
t211 = t105.*5.902360448889528e-2;
t212 = t106.*5.902360448889528e-2;
t213 = t105.*7.043523830355962e-1;
t214 = t106.*7.043523830355962e-1;
t216 = t111.*6.966785627840228e-2;
t217 = t110.*6.966785627840228e-2;
t218 = t111.*2.078572828445888e-1;
t219 = t110.*2.078572828445888e-1;
t220 = t111.*7.17599677487021e-1;
t221 = t110.*7.17599677487021e-1;
t222 = t111.*1.096343788571921e-1;
t223 = t110.*1.096343788571921e-1;
t224 = t73.*1.353803736954536;
t225 = t72.*1.353803736954536;
t226 = t73.*6.769018684772682e-1;
t227 = t72.*1.06229898470213e+1;
t228 = t73.*1.06229898470213e+1;
t229 = t77.*2.857199915598597;
t230 = t77.*5.714399831197195;
t231 = t72.*6.769018684772682e-1;
t232 = t78.*2.857199915598597;
t233 = t78.*5.714399831197195;
t234 = t72.*4.241794227085137;
t235 = t73.*4.241794227085137;
t236 = t72.*5.311494923510651;
t237 = t73.*5.311494923510651;
t239 = t72.*1.184934912206368;
t240 = t72.*2.369869824412735;
t241 = t73.*1.184934912206368;
t242 = t73.*2.369869824412735;
t244 = t72.*2.120897113542568;
t246 = t73.*2.120897113542568;
t250 = t77.*1.294704791288842;
t251 = t77.*2.589409582577684;
t252 = t78.*1.294704791288842;
t253 = t78.*2.589409582577684;
t254 = t77.*1.147826242706024;
t255 = t78.*1.147826242706024;
t259 = t77.*5.739131213530122e-1;
t260 = t77.*2.001174397900496e+1;
t261 = t77.*4.002348795800993e+1;
t263 = t78.*5.739131213530122e-1;
t264 = t78.*2.001174397900496e+1;
t265 = t78.*4.002348795800993e+1;
t266 = t84.*1.968054798479955;
t267 = t84.*5.348406989475424;
t268 = t85.*1.968054798479955;
t269 = t85.*5.348406989475424;
t270 = t95.*3.306284947566274;
t271 = t95.*6.612569895132548;
t272 = t89.*3.306284947566274;
t273 = t89.*6.612569895132548;
t275 = t84.*9.840273992399776e-1;
t277 = t84.*2.674203494737712;
t278 = t85.*9.840273992399776e-1;
t280 = t85.*2.674203494737712;
t281 = t84.*4.751064619307353e+1;
t282 = t85.*4.751064619307353e+1;
t285 = t95.*1.39342291567647e+1;
t286 = t95.*2.786845831352939e+1;
t287 = t89.*1.39342291567647e+1;
t288 = t89.*2.786845831352939e+1;
t290 = t84.*3.521663908435019;
t291 = t84.*7.043327816870038;
t294 = t85.*3.521663908435019;
t295 = t85.*7.043327816870038;
t296 = t95.*1.894171845065846;
t297 = t95.*3.788343690131691;
t299 = t89.*1.894171845065846;
t300 = t89.*3.788343690131691;
t301 = t84.*2.375532309653677e+1;
t302 = t85.*2.375532309653677e+1;
t303 = t89.*1.621251186419005;
t304 = t89.*3.24250237283801;
t305 = t95.*1.621251186419005;
t306 = t95.*3.24250237283801;
t309 = t105.*1.243239735437342e+1;
t310 = t106.*1.243239735437342e+1;
t311 = t105.*2.141972367605808;
t312 = t105.*4.283944735211616;
t313 = t106.*2.141972367605808;
t314 = t106.*4.283944735211616;
t315 = t105.*6.216198677186712;
t316 = t106.*6.216198677186712;
t317 = t105.*2.556101672515437e+1;
t318 = t105.*5.112203345030875e+1;
t319 = t105.*3.876637815316351;
t320 = t105.*7.753275630632703;
t321 = t106.*2.556101672515437e+1;
t322 = t106.*5.112203345030875e+1;
t323 = t106.*3.876637815316351;
t324 = t106.*7.753275630632703;
t328 = t111.*2.528253303925992;
t329 = t110.*2.528253303925992;
t330 = t111.*7.543161081300932;
t331 = t110.*7.543161081300932;
t332 = t111.*2.604176233373252e+1;
t333 = t110.*2.604176233373252e+1;
t334 = t111.*3.97864230904288;
t335 = t110.*3.97864230904288;
t354 = t102.*6.491387775649387e-1;
t362 = t102.*4.875659379893039;
t369 = t102.*8.790578612396049e-1;
t378 = t102.*8.289508866590401e-1;
t379 = t102.*1.494558451938363e-1;
t383 = t102.*1.12255779300794;
t390 = t72.*t76.*7.386398980558217e-3;
t391 = t73.*t76.*7.386398980558217e-3;
t392 = t72.*t75.*8.210942169983652e-2;
t393 = t73.*t75.*8.210942169983652e-2;
t394 = t72.*t76.*5.795939192266646e-2;
t395 = t73.*t76.*5.795939192266646e-2;
t396 = t72.*t75.*4.154570279669746e-2;
t397 = t73.*t75.*4.154570279669746e-2;
t400 = t72.*t75.*5.55465297295666e-2;
t401 = t73.*t75.*5.55465297295666e-2;
t402 = t72.*t75.*6.20217245890067e-2;
t403 = t73.*t75.*6.20217245890067e-2;
t404 = t72.*t76.*1.293008992165735e-2;
t405 = t73.*t76.*1.293008992165735e-2;
t406 = t72.*t76.*2.314337277954461e-2;
t407 = t73.*t76.*2.314337277954461e-2;
t410 = t77.*t101.*3.117795876572799e-2;
t411 = t78.*t101.*3.117795876572799e-2;
t412 = t77.*t100.*7.873932755209373e-2;
t413 = t87.*t95.*1.810676324582379e-1;
t414 = t78.*t100.*7.873932755209373e-2;
t415 = t87.*t89.*1.810676324582379e-1;
t421 = t77.*t101.*6.262578804852617e-3;
t422 = t78.*t101.*6.262578804852617e-3;
t423 = t77.*t101.*1.412790626802782e-2;
t424 = t78.*t101.*1.412790626802782e-2;
t425 = t77.*t101.*2.183695040733637e-1;
t426 = t78.*t101.*2.183695040733637e-1;
t427 = t84.*t98.*3.700264708997962e-2;
t428 = t85.*t98.*3.700264708997962e-2;
t429 = t78.*t100.*1.993938186310297e-1;
t431 = t88.*t95.*3.607840501471245e-2;
t432 = t88.*t89.*3.607840501471245e-2;
t435 = t87.*t95.*1.348317732682517e-1;
t436 = t77.*t100.*2.859465635352065e-1;
t437 = t87.*t89.*1.348317732682517e-1;
t438 = t77.*t100.*1.993938186310297e-1;
t439 = t78.*t100.*2.859465635352065e-1;
t440 = t88.*t95.*2.0669331312189e-2;
t441 = t88.*t89.*2.0669331312189e-2;
t442 = t84.*t99.*2.91810894417692e-2;
t444 = t85.*t99.*2.91810894417692e-2;
t445 = t88.*t95.*1.520512511952798e-1;
t446 = t84.*t99.*1.073777354897583e-2;
t447 = t88.*t89.*1.520512511952798e-1;
t448 = t85.*t99.*1.073777354897583e-2;
t449 = t87.*t95.*1.573732608437441e-1;
t450 = t88.*t95.*1.769120262222513e-2;
t451 = t87.*t89.*1.573732608437441e-1;
t452 = t84.*t98.*3.694071270429475e-1;
t453 = t88.*t89.*1.769120262222513e-2;
t454 = t85.*t98.*3.694071270429475e-1;
t455 = t87.*t95.*1.231273430743841e-1;
t456 = t87.*t89.*1.231273430743841e-1;
t460 = t77.*t100.*1.345824284957329e-1;
t461 = t84.*t99.*3.842863480588345e-2;
t462 = t84.*t98.*5.035148185091873e-2;
t464 = t78.*t100.*1.345824284957329e-1;
t465 = t85.*t99.*3.842863480588345e-2;
t466 = t85.*t98.*5.035148185091873e-2;
t474 = t84.*t99.*2.592196926532531e-1;
t475 = t85.*t99.*2.592196926532531e-1;
t476 = t84.*t98.*2.6526615661289e-1;
t477 = t85.*t98.*2.6526615661289e-1;
t480 = t92.*t106.*3.149797662150347e-1;
t481 = t92.*t105.*5.678933207774932e-2;
t482 = t92.*t106.*5.678933207774932e-2;
t483 = t92.*t105.*3.149797662150347e-1;
t484 = t93.*t105.*2.337334737760253e-2;
t485 = t93.*t106.*2.337334737760253e-2;
t486 = t90.*t111.*8.907876918584193e-4;
t487 = t93.*t105.*4.230213409140225e-2;
t488 = t90.*t110.*8.907876918584193e-4;
t489 = t93.*t106.*4.230213409140225e-2;
t490 = t93.*t105.*6.783158048508356e-2;
t491 = t93.*t106.*6.783158048508356e-2;
t492 = t92.*t105.*4.26542750475325e-1;
t493 = t92.*t106.*4.26542750475325e-1;
t494 = t93.*t105.*2.789235436820961e-1;
t495 = t93.*t106.*2.789235436820961e-1;
t501 = t92.*t105.*7.252003547718363e-2;
t502 = t92.*t106.*7.252003547718363e-2;
t503 = t91.*t111.*2.75884710862473e-2;
t504 = t91.*t110.*2.75884710862473e-2;
t507 = t91.*t111.*8.231148400645715e-2;
t508 = t91.*t110.*8.231148400645715e-2;
t509 = t90.*t111.*3.426310303908123e-1;
t510 = t90.*t110.*3.426310303908123e-1;
t511 = t91.*t111.*2.841694722848603e-1;
t513 = t91.*t110.*2.841694722848603e-1;
t514 = t90.*t111.*4.550026298466454e-1;
t515 = t90.*t110.*4.550026298466454e-1;
t516 = t90.*t111.*1.626186884310744e-1;
t517 = t90.*t110.*1.626186884310744e-1;
t519 = t91.*t111.*4.341521402744809e-2;
t520 = t91.*t110.*4.341521402744809e-2;
t524 = t72.*t76.*2.103351989710218;
t525 = t73.*t76.*2.103351989710218;
t526 = t72.*t76.*2.680531399169982e-1;
t527 = t73.*t76.*2.680531399169982e-1;
t528 = t72.*t75.*2.979758927366625;
t529 = t73.*t75.*2.979758927366625;
t530 = t72.*t75.*1.507697609352744;
t531 = t73.*t75.*1.507697609352744;
t534 = t72.*t76.*8.398752569628571e-1;
t535 = t72.*t75.*2.250774438655373;
t536 = t73.*t76.*8.398752569628571e-1;
t537 = t72.*t75.*2.015788985227274;
t538 = t73.*t75.*2.250774438655373;
t539 = t73.*t75.*2.015788985227274;
t540 = t72.*t76.*4.692342252337216e-1;
t541 = t73.*t76.*4.692342252337216e-1;
t544 = t77.*t101.*1.131451166577045;
t545 = t78.*t101.*1.131451166577045;
t547 = t77.*t100.*2.857457881823851;
t548 = t78.*t100.*2.857457881823851;
t549 = t87.*t95.*6.570962054110382;
t550 = t87.*t89.*6.570962054110382;
t553 = t77.*t101.*7.924650615685965;
t554 = t77.*t101.*2.272695960557928e-1;
t556 = t78.*t101.*7.924650615685965;
t557 = t78.*t101.*2.272695960557928e-1;
t559 = t77.*t101.*5.127030973503815e-1;
t560 = t78.*t101.*5.127030973503815e-1;
t563 = t84.*t98.*1.342829674353717;
t564 = t87.*t95.*4.893058211485927;
t565 = t85.*t98.*1.342829674353717;
t566 = t87.*t89.*4.893058211485927;
t568 = t88.*t89.*1.309288839236244;
t569 = t88.*t95.*7.500920506460749e-1;
t570 = t88.*t89.*7.500920506460749e-1;
t571 = t77.*t100.*7.236021138956767;
t572 = t78.*t100.*7.236021138956767;
t573 = t77.*t100.*1.037702869907725e+1;
t574 = t88.*t95.*1.309288839236244;
t575 = t78.*t100.*1.037702869907725e+1;
t578 = t84.*t98.*1.340582069452416e+1;
t579 = t85.*t98.*1.340582069452416e+1;
t580 = t84.*t99.*1.058984583916134;
t581 = t85.*t99.*1.058984583916134;
t582 = t88.*t95.*5.51795474607882;
t583 = t84.*t99.*3.896748500990311e-1;
t584 = t88.*t89.*5.51795474607882;
t585 = t85.*t99.*3.896748500990311e-1;
t586 = t87.*t95.*4.468303297398084;
t588 = t87.*t95.*5.711090995649731;
t589 = t87.*t89.*4.468303297398084;
t591 = t87.*t89.*5.711090995649731;
t592 = t84.*t99.*9.407107946228559;
t593 = t85.*t99.*1.394578907740267;
t594 = t85.*t99.*9.407107946228559;
t597 = t88.*t95.*6.420154698219261e-1;
t598 = t88.*t89.*6.420154698219261e-1;
t600 = t77.*t100.*4.884009465355168;
t601 = t78.*t100.*4.884009465355168;
t604 = t84.*t98.*9.626534713458666;
t605 = t85.*t98.*9.626534713458666;
t607 = t84.*t99.*1.394578907740267;
t612 = t84.*t98.*1.82726019067447;
t613 = t85.*t98.*1.82726019067447;
t614 = t92.*t105.*2.060890403740334;
t615 = t92.*t106.*2.060890403740334;
t616 = t92.*t105.*1.143064645796879e+1;
t617 = t92.*t106.*1.143064645796879e+1;
t618 = t93.*t105.*8.482210575719e-1;
t619 = t93.*t106.*8.482210575719e-1;
t620 = t92.*t105.*1.547927804532199e+1;
t621 = t92.*t106.*1.547927804532199e+1;
t622 = t93.*t105.*2.461614676165938;
t623 = t93.*t106.*2.461614676165938;
t624 = t90.*t111.*3.232677227842077e-2;
t625 = t90.*t110.*3.232677227842077e-2;
t629 = t93.*t105.*1.012216262316113e+1;
t630 = t93.*t105.*1.535148574865275;
t631 = t93.*t106.*1.012216262316113e+1;
t632 = t93.*t106.*1.535148574865275;
t633 = t92.*t105.*2.631759165422457;
t634 = t92.*t106.*2.631759165422457;
t637 = t91.*t111.*1.001188308354693;
t638 = t91.*t110.*1.001188308354693;
t641 = t90.*t111.*1.243411353367114e+1;
t642 = t90.*t110.*1.243411353367114e+1;
t643 = t90.*t111.*1.651208984539144e+1;
t644 = t90.*t110.*1.651208984539144e+1;
t645 = t91.*t111.*2.987091788195169;
t646 = t91.*t110.*2.987091788195169;
t647 = t91.*t110.*1.031253788415808e+1;
t651 = t91.*t111.*1.031253788415808e+1;
t653 = t90.*t111.*5.901448074747681;
t654 = t90.*t110.*5.901448074747681;
t655 = t91.*t111.*1.57554235438098;
t656 = t91.*t110.*1.57554235438098;
t126 = -t120;
t131 = -t125;
t132 = -t124;
t133 = -t128;
t134 = -t129;
t136 = -t130;
t146 = -t139;
t147 = -t141;
t148 = -t143;
t149 = -t145;
t158 = -t155;
t159 = -t154;
MMred = ft_1({t102,t104,t107,t108,t112,t113,t114,t115,t116,t117,t118,t119,t121,t122,t123,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149,t150,t151,t152,t153,t155,t156,t157,t158,t159,t161,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t173,t174,t175,t178,t183,t184,t185,t186,t187,t188,t189,t190,t191,t192,t193,t194,t197,t198,t199,t200,t201,t202,t204,t206,t207,t208,t209,t210,t211,t212,t213,t214,t216,t217,t218,t219,t220,t221,t222,t223,t224,t225,t226,t227,t228,t229,t230,t231,t232,t233,t234,t235,t236,t237,t239,t240,t241,t242,t244,t246,t250,t251,t252,t253,t254,t255,t259,t260,t261,t263,t264,t265,t266,t267,t268,t269,t270,t271,t272,t273,t275,t277,t278,t280,t281,t282,t285,t286,t287,t288,t290,t291,t294,t295,t296,t297,t299,t300,t301,t302,t303,t304,t305,t306,t309,t310,t311,t312,t313,t314,t315,t316,t317,t318,t319,t320,t321,t322,t323,t324,t328,t329,t330,t331,t332,t333,t334,t335,t354,t362,t369,t378,t379,t383,t390,t391,t392,t393,t394,t395,t396,t397,t400,t401,t402,t403,t404,t405,t406,t407,t410,t411,t412,t413,t414,t415,t421,t422,t423,t424,t425,t426,t427,t428,t429,t431,t432,t435,t436,t437,t438,t439,t440,t441,t442,t444,t445,t446,t447,t448,t449,t450,t451,t452,t453,t454,t455,t456,t460,t461,t462,t464,t465,t466,t474,t475,t476,t477,t480,t481,t482,t483,t484,t485,t486,t487,t488,t489,t490,t491,t492,t493,t494,t495,t501,t502,t503,t504,t507,t508,t509,t510,t511,t513,t514,t515,t516,t517,t519,t520,t524,t525,t526,t527,t528,t529,t530,t531,t534,t535,t536,t537,t538,t539,t540,t541,t544,t545,t547,t548,t549,t550,t553,t554,t556,t557,t559,t560,t563,t564,t565,t566,t568,t569,t570,t571,t572,t573,t574,t575,t578,t579,t580,t581,t582,t583,t584,t585,t586,t588,t589,t591,t592,t593,t594,t597,t598,t600,t601,t604,t605,t607,t612,t613,t614,t615,t616,t617,t618,t619,t620,t621,t622,t623,t624,t625,t629,t630,t631,t632,t633,t634,t637,t638,t641,t642,t643,t644,t645,t646,t647,t651,t653,t654,t655,t656,t79,t97});
end
function MMred = ft_1(ct)
t102 = ct{1};
t104 = ct{2};
t107 = ct{3};
t108 = ct{4};
t112 = ct{5};
t113 = ct{6};
t114 = ct{7};
t115 = ct{8};
t116 = ct{9};
t117 = ct{10};
t118 = ct{11};
t119 = ct{12};
t121 = ct{13};
t122 = ct{14};
t123 = ct{15};
t125 = ct{16};
t126 = ct{17};
t127 = ct{18};
t128 = ct{19};
t129 = ct{20};
t130 = ct{21};
t131 = ct{22};
t132 = ct{23};
t133 = ct{24};
t134 = ct{25};
t135 = ct{26};
t136 = ct{27};
t137 = ct{28};
t138 = ct{29};
t139 = ct{30};
t140 = ct{31};
t141 = ct{32};
t142 = ct{33};
t143 = ct{34};
t144 = ct{35};
t145 = ct{36};
t146 = ct{37};
t147 = ct{38};
t148 = ct{39};
t149 = ct{40};
t150 = ct{41};
t151 = ct{42};
t152 = ct{43};
t153 = ct{44};
t155 = ct{45};
t156 = ct{46};
t157 = ct{47};
t158 = ct{48};
t159 = ct{49};
t161 = ct{50};
t162 = ct{51};
t163 = ct{52};
t164 = ct{53};
t165 = ct{54};
t166 = ct{55};
t167 = ct{56};
t168 = ct{57};
t169 = ct{58};
t170 = ct{59};
t171 = ct{60};
t172 = ct{61};
t173 = ct{62};
t174 = ct{63};
t175 = ct{64};
t178 = ct{65};
t183 = ct{66};
t184 = ct{67};
t185 = ct{68};
t186 = ct{69};
t187 = ct{70};
t188 = ct{71};
t189 = ct{72};
t190 = ct{73};
t191 = ct{74};
t192 = ct{75};
t193 = ct{76};
t194 = ct{77};
t197 = ct{78};
t198 = ct{79};
t199 = ct{80};
t200 = ct{81};
t201 = ct{82};
t202 = ct{83};
t204 = ct{84};
t206 = ct{85};
t207 = ct{86};
t208 = ct{87};
t209 = ct{88};
t210 = ct{89};
t211 = ct{90};
t212 = ct{91};
t213 = ct{92};
t214 = ct{93};
t216 = ct{94};
t217 = ct{95};
t218 = ct{96};
t219 = ct{97};
t220 = ct{98};
t221 = ct{99};
t222 = ct{100};
t223 = ct{101};
t224 = ct{102};
t225 = ct{103};
t226 = ct{104};
t227 = ct{105};
t228 = ct{106};
t229 = ct{107};
t230 = ct{108};
t231 = ct{109};
t232 = ct{110};
t233 = ct{111};
t234 = ct{112};
t235 = ct{113};
t236 = ct{114};
t237 = ct{115};
t239 = ct{116};
t240 = ct{117};
t241 = ct{118};
t242 = ct{119};
t244 = ct{120};
t246 = ct{121};
t250 = ct{122};
t251 = ct{123};
t252 = ct{124};
t253 = ct{125};
t254 = ct{126};
t255 = ct{127};
t259 = ct{128};
t260 = ct{129};
t261 = ct{130};
t263 = ct{131};
t264 = ct{132};
t265 = ct{133};
t266 = ct{134};
t267 = ct{135};
t268 = ct{136};
t269 = ct{137};
t270 = ct{138};
t271 = ct{139};
t272 = ct{140};
t273 = ct{141};
t275 = ct{142};
t277 = ct{143};
t278 = ct{144};
t280 = ct{145};
t281 = ct{146};
t282 = ct{147};
t285 = ct{148};
t286 = ct{149};
t287 = ct{150};
t288 = ct{151};
t290 = ct{152};
t291 = ct{153};
t294 = ct{154};
t295 = ct{155};
t296 = ct{156};
t297 = ct{157};
t299 = ct{158};
t300 = ct{159};
t301 = ct{160};
t302 = ct{161};
t303 = ct{162};
t304 = ct{163};
t305 = ct{164};
t306 = ct{165};
t309 = ct{166};
t310 = ct{167};
t311 = ct{168};
t312 = ct{169};
t313 = ct{170};
t314 = ct{171};
t315 = ct{172};
t316 = ct{173};
t317 = ct{174};
t318 = ct{175};
t319 = ct{176};
t320 = ct{177};
t321 = ct{178};
t322 = ct{179};
t323 = ct{180};
t324 = ct{181};
t328 = ct{182};
t329 = ct{183};
t330 = ct{184};
t331 = ct{185};
t332 = ct{186};
t333 = ct{187};
t334 = ct{188};
t335 = ct{189};
t354 = ct{190};
t362 = ct{191};
t369 = ct{192};
t378 = ct{193};
t379 = ct{194};
t383 = ct{195};
t390 = ct{196};
t391 = ct{197};
t392 = ct{198};
t393 = ct{199};
t394 = ct{200};
t395 = ct{201};
t396 = ct{202};
t397 = ct{203};
t400 = ct{204};
t401 = ct{205};
t402 = ct{206};
t403 = ct{207};
t404 = ct{208};
t405 = ct{209};
t406 = ct{210};
t407 = ct{211};
t410 = ct{212};
t411 = ct{213};
t412 = ct{214};
t413 = ct{215};
t414 = ct{216};
t415 = ct{217};
t421 = ct{218};
t422 = ct{219};
t423 = ct{220};
t424 = ct{221};
t425 = ct{222};
t426 = ct{223};
t427 = ct{224};
t428 = ct{225};
t429 = ct{226};
t431 = ct{227};
t432 = ct{228};
t435 = ct{229};
t436 = ct{230};
t437 = ct{231};
t438 = ct{232};
t439 = ct{233};
t440 = ct{234};
t441 = ct{235};
t442 = ct{236};
t444 = ct{237};
t445 = ct{238};
t446 = ct{239};
t447 = ct{240};
t448 = ct{241};
t449 = ct{242};
t450 = ct{243};
t451 = ct{244};
t452 = ct{245};
t453 = ct{246};
t454 = ct{247};
t455 = ct{248};
t456 = ct{249};
t460 = ct{250};
t461 = ct{251};
t462 = ct{252};
t464 = ct{253};
t465 = ct{254};
t466 = ct{255};
t474 = ct{256};
t475 = ct{257};
t476 = ct{258};
t477 = ct{259};
t480 = ct{260};
t481 = ct{261};
t482 = ct{262};
t483 = ct{263};
t484 = ct{264};
t485 = ct{265};
t486 = ct{266};
t487 = ct{267};
t488 = ct{268};
t489 = ct{269};
t490 = ct{270};
t491 = ct{271};
t492 = ct{272};
t493 = ct{273};
t494 = ct{274};
t495 = ct{275};
t501 = ct{276};
t502 = ct{277};
t503 = ct{278};
t504 = ct{279};
t507 = ct{280};
t508 = ct{281};
t509 = ct{282};
t510 = ct{283};
t511 = ct{284};
t513 = ct{285};
t514 = ct{286};
t515 = ct{287};
t516 = ct{288};
t517 = ct{289};
t519 = ct{290};
t520 = ct{291};
t524 = ct{292};
t525 = ct{293};
t526 = ct{294};
t527 = ct{295};
t528 = ct{296};
t529 = ct{297};
t530 = ct{298};
t531 = ct{299};
t534 = ct{300};
t535 = ct{301};
t536 = ct{302};
t537 = ct{303};
t538 = ct{304};
t539 = ct{305};
t540 = ct{306};
t541 = ct{307};
t544 = ct{308};
t545 = ct{309};
t547 = ct{310};
t548 = ct{311};
t549 = ct{312};
t550 = ct{313};
t553 = ct{314};
t554 = ct{315};
t556 = ct{316};
t557 = ct{317};
t559 = ct{318};
t560 = ct{319};
t563 = ct{320};
t564 = ct{321};
t565 = ct{322};
t566 = ct{323};
t568 = ct{324};
t569 = ct{325};
t570 = ct{326};
t571 = ct{327};
t572 = ct{328};
t573 = ct{329};
t574 = ct{330};
t575 = ct{331};
t578 = ct{332};
t579 = ct{333};
t580 = ct{334};
t581 = ct{335};
t582 = ct{336};
t583 = ct{337};
t584 = ct{338};
t585 = ct{339};
t586 = ct{340};
t588 = ct{341};
t589 = ct{342};
t591 = ct{343};
t592 = ct{344};
t593 = ct{345};
t594 = ct{346};
t597 = ct{347};
t598 = ct{348};
t600 = ct{349};
t601 = ct{350};
t604 = ct{351};
t605 = ct{352};
t607 = ct{353};
t612 = ct{354};
t613 = ct{355};
t614 = ct{356};
t615 = ct{357};
t616 = ct{358};
t617 = ct{359};
t618 = ct{360};
t619 = ct{361};
t620 = ct{362};
t621 = ct{363};
t622 = ct{364};
t623 = ct{365};
t624 = ct{366};
t625 = ct{367};
t629 = ct{368};
t630 = ct{369};
t631 = ct{370};
t632 = ct{371};
t633 = ct{372};
t634 = ct{373};
t637 = ct{374};
t638 = ct{375};
t641 = ct{376};
t642 = ct{377};
t643 = ct{378};
t644 = ct{379};
t645 = ct{380};
t646 = ct{381};
t647 = ct{382};
t651 = ct{383};
t653 = ct{384};
t654 = ct{385};
t655 = ct{386};
t656 = ct{387};
t79 = ct{388};
t97 = ct{389};
t160 = -t156;
t176 = -t163;
t177 = -t162;
t179 = -t166;
t180 = -t165;
t181 = -t168;
t182 = -t170;
t195 = -t190;
t196 = -t194;
t203 = -t197;
t205 = -t199;
t215 = -t211;
t238 = -t231;
t243 = -t234;
t245 = -t235;
t247 = -t240;
t248 = -t242;
t249 = -t244;
t256 = -t251;
t257 = -t253;
t258 = -t254;
t262 = -t255;
t274 = -t266;
t276 = -t267;
t279 = -t269;
t283 = -t271;
t284 = -t273;
t289 = -t275;
t292 = -t277;
t293 = -t278;
t298 = -t280;
t307 = -t304;
t308 = -t306;
t325 = -t311;
t326 = -t315;
t327 = -t316;
t336 = t79.*1.655152678459117e-1;
t337 = t79.*1.237961788021798e-1;
t338 = t79.*1.848097875347668e-1;
t339 = t79.*8.374737008110868e-2;
t340 = t97.*8.859744658340967e-1;
t341 = t79.*9.351000589077968e-2;
t346 = t79.*1.250227092713193e-1;
t347 = t97.*8.090650991652088e-1;
t348 = t97.*1.034093725335858;
t351 = t97.*7.700365262397986e-1;
t352 = t104.*1.107615529496445e-2;
t353 = t108.*8.170802617541204e-1;
t355 = t108.*5.697594386341912e-1;
t356 = t107.*4.960508521110232e-1;
t357 = t104.*1.470876640113508e-2;
t360 = t97.*6.0246925708856e-1;
t363 = t107.*3.562072667015491e-1;
t364 = t97.*7.031914603503905e-1;
t366 = t104.*5.256937309126799e-3;
t368 = t107.*6.761346397709627e-2;
t370 = t108.*3.845636210601268e-1;
t374 = t107.*3.556110531877969;
t375 = t104.*5.657554357632136;
t376 = t108.*2.069115378352799;
t377 = t107.*6.750029373970042e-1;
t381 = t107.*4.847102879119532e-1;
t382 = t104.*2.022019234648673;
t385 = t104.*2.685174393913617;
t386 = t108.*1.396565722891732;
t387 = t108.*9.738412975272607e-1;
t398 = -t393;
t399 = -t395;
t408 = -t403;
t409 = -t405;
t416 = -t411;
t417 = -t412;
t418 = -t413;
t419 = -t414;
t420 = -t415;
t430 = -t422;
t433 = -t424;
t434 = -t426;
t443 = -t432;
t457 = -t436;
t458 = -t439;
t459 = -t441;
t463 = -t442;
t467 = -t446;
t468 = -t447;
t469 = -t449;
t470 = -t451;
t471 = -t452;
t472 = -t453;
t473 = -t454;
t478 = -t465;
t479 = -t475;
t496 = -t489;
t497 = -t490;
t498 = -t492;
t499 = -t493;
t500 = -t495;
t505 = -t501;
t506 = -t504;
t512 = -t508;
t518 = -t513;
t521 = -t514;
t522 = -t515;
t523 = -t520;
t532 = -t525;
t533 = -t529;
t542 = -t538;
t543 = -t541;
t546 = -t545;
t551 = -t547;
t552 = -t548;
t555 = -t549;
t558 = -t550;
t561 = -t556;
t562 = -t557;
t567 = -t560;
t576 = -t568;
t577 = -t570;
t587 = -t573;
t590 = -t575;
t595 = -t578;
t596 = -t579;
t599 = -t580;
t602 = -t583;
t603 = -t584;
t606 = -t588;
t608 = -t591;
t609 = -t593;
t610 = -t594;
t611 = -t598;
t626 = -t620;
t627 = -t621;
t628 = -t622;
t635 = -t631;
t636 = -t632;
t639 = -t633;
t640 = -t638;
t648 = -t643;
t649 = -t644;
t650 = -t646;
t652 = -t647;
t657 = -t656;
t658 = t123+t394+t397;
t664 = t137+t401+t404;
t666 = t237+t524+t531;
t672 = t241+t539+t540;
t676 = t117+t185+t445+t456;
t680 = t130+t189+t437+t450;
t685 = t228+t287+t582+t589;
t688 = t242+t303+t566+t597;
t694 = t117+t152+t186+t425+t429;
t696 = t130+t144+t190+t421+t464;
t702 = t228+t264+t288+t553+t572;
t704 = t242+t263+t304+t554+t601;
t710 = t117+t153+t173+t186+t474+t477;
t718 = t228+t265+t288+t302+t592+t605;
t726 = t117+t153+t174+t186+t214+t480+t494;
t738 = t228+t265+t282+t288+t321+t617+t629;
t739 = t117+t153+t174+t186+t206+t221+t510+t511;
t748 = t228+t265+t282+t288+t322+t333+t642+t651;
t342 = -t337;
t343 = -t338;
t344 = -t339;
t345 = -t340;
t349 = -t346;
t350 = -t348;
t358 = -t352;
t359 = -t353;
t361 = -t355;
t365 = -t357;
t367 = -t360;
t371 = -t364;
t372 = -t366;
t373 = -t368;
t380 = -t377;
t384 = -t381;
t388 = -t386;
t389 = -t387;
t659 = t126+t391+t402;
t660 = t121+t390+t408;
t661 = t122+t396+t399;
t662 = t132+t392+t407;
t663 = t127+t398+t406;
t665 = t135+t400+t409;
t667 = t236+t530+t532;
t668 = t238+t527+t535;
t669 = t226+t526+t542;
t670 = t249+t528+t536;
t671 = t246+t533+t534;
t673 = t239+t537+t543;
t674 = t115+t192+t420+t440;
t675 = t114+t191+t418+t459;
t677 = t128+t169+t431+t470;
t678 = t116+t183+t455+t468;
t679 = t125+t167+t443+t469;
t681 = t129+t193+t435+t472;
t682 = t224+t299+t558+t569;
t683 = t225+t296+t555+t577;
t684 = t235+t272+t574+t608;
t686 = t227+t285+t586+t603;
t687 = t234+t270+t576+t606;
t689 = t240+t305+t564+t611;
t690 = t115+t119+t188+t410+t458;
t691 = t114+t118+t187+t416+t457;
t692 = t128+t140+t170+t419+t423;
t693 = t125+t138+t168+t417+t433;
t695 = t116+t150+t184+t434+t438;
t697 = t129+t142+t194+t430+t460;
t698 = t224+t232+t300+t544+t590;
t699 = t225+t229+t297+t546+t587;
t700 = t235+t252+t273+t552+t559;
t701 = t234+t250+t271+t551+t567;
t703 = t227+t260+t286+t561+t571;
t705 = t240+t259+t306+t562+t600;
t706 = t113+t115+t178+t188+t461+t473;
t707 = t125+t139+t168+t177+t427+t444;
t708 = t112+t114+t175+t187+t471+t478;
t709 = t128+t141+t170+t180+t428+t463;
t711 = t129+t143+t159+t194+t448+t462;
t712 = t116+t151+t171+t184+t476+t479;
t713 = t130+t145+t160+t190+t466+t467;
t714 = t224+t233+t294+t300+t596+t607;
t715 = t234+t251+t271+t292+t563+t581;
t716 = t235+t253+t273+t298+t565+t599;
t717 = t225+t230+t290+t297+t595+t609;
t719 = t227+t261+t286+t301+t604+t610;
t720 = t240+t254+t289+t306+t585+t612;
t721 = t242+t255+t293+t304+t602+t613;
t722 = t113+t115+t164+t188+t209+t487+t499;
t723 = t112+t114+t161+t187+t207+t496+t498;
t724 = t125+t139+t168+t176+t203+t481+t491;
t725 = t128+t141+t170+t179+t205+t482+t497;
t727 = t116+t151+t172+t184+t213+t483+t500;
t728 = t129+t143+t158+t194+t215+t485+t505;
t729 = t136+t149+t157+t195+t212+t484+t502;
t730 = t133+t147+t166+t182+t200+t219+t488+t507;
t731 = t131+t146+t163+t181+t198+t218+t486+t512;
t732 = t224+t233+t295+t300+t323+t627+t630;
t733 = t113+t115+t164+t188+t210+t223+t519+t522;
t734 = t234+t251+t271+t276+t326+t614+t623;
t735 = t112+t114+t161+t187+t208+t222+t521+t523;
t736 = t225+t230+t291+t297+t319+t626+t636;
t737 = t136+t149+t157+t195+t202+t217+t503+t517;
t740 = t235+t253+t273+t279+t327+t615+t628;
t741 = t134+t148+t155+t196+t201+t216+t506+t516;
t742 = t227+t261+t281+t286+t317+t616+t635;
t743 = t116+t151+t172+t184+t204+t220+t509+t518;
t744 = t240+t254+t274+t306+t325+t619+t639;
t745 = t248+t262+t268+t307+t313+t618+t634;
t746 = t245+t257+t269+t284+t310+t331+t625+t645;
t747 = t243+t256+t267+t283+t309+t330+t624+t650;
t749 = t227+t261+t281+t286+t318+t332+t641+t652;
t750 = t224+t233+t295+t300+t324+t335+t649+t655;
t751 = t225+t230+t291+t297+t320+t334+t648+t657;
t752 = t248+t262+t268+t307+t314+t329+t637+t654;
t753 = t247+t258+t266+t308+t312+t328+t640+t653;
et1 = t79.*6.263835678546036e-2;
et2 = t97.*5.501703130591357e-1;
et3 = t102.*3.600422349017724;
et4 = t104.*4.260313132038088;
et5 = t107.*2.553593864939749;
et6 = t108.*1.442817886591442+t658.*t666+t661.*t667+t676.*t685+t678.*t686+t694.*t702+t695.*t703+t710.*t718+t712.*t719+t726.*t738+t727.*t742+t739.*t748+t743.*t749;
et7 = t79.*1.395969123462597e-1;
et8 = t97.*1.189788542110687;
et9 = t102.*6.602573832823954;
et10 = t104.*7.513044303917191;
et11 = t107.*4.952205708417063;
et12 = t108.*2.967275696207353+t659.*t668+t660.*t669+t674.*t682+t675.*t683+t690.*t698+t691.*t699+t706.*t714+t708.*t717+t722.*t732+t723.*t736+t733.*t750+t735.*t751;
et13 = t79.*1.119700827964597e-1;
et14 = t97.*6.597397153594278e-1;
et15 = t102.*1.908552680438397e-1;
et16 = t104.*9.596857457595569e-1;
et17 = t107.*9.200525832765186e-2;
et18 = t108.*6.573018546436445e-1+(-t129+t148+t155+t196+t201+t216+t506+t516).*(-t240+t258+t266+t308+t312+t328+t640+t653)+t664.*t672+t665.*t673+t680.*t688+t681.*t689+t696.*t704+t697.*t705+t711.*t720+t713.*t721+t728.*t744+t729.*t745+t737.*t752;
et19 = t79.*2.446662823309988e-1;
et20 = t97.*8.987730129607432e-1;
et21 = t102.*1.170365895138567e-1;
et22 = t104.*2.879629086312717e-5;
et23 = t107.*4.968825254206283e-2;
et24 = t108.*2.249943121232401e-1+(-t125+t146+t163+t181+t198+t218+t486+t512).*(-t234+t256+t267+t283+t309+t330+t624+t650)+t662.*t670+t663.*t671+t677.*t684+t679.*t687+t692.*t700+t693.*t701+t707.*t715+t709.*t716+t724.*t734+t725.*t740+t730.*t746;
mt1 = [et1+et2+et3+et4+et5+et6,t341+t347+t362+t374+t375+t376-t658.*t669+t661.*t668-t676.*t682-t678.*t683-t694.*t698-t695.*t699-t710.*t714-t712.*t717-t726.*t732-t727.*t736-t739.*t750-t743.*t751,t344+t367+t378+t382+t384+t389+t743.*(-t240+t258+t266+t308+t312+t328+t640+t653)-t658.*t672-t661.*t673-t676.*t688-t678.*t689-t694.*t704-t695.*t705-t710.*t721-t712.*t720+t726.*t745-t727.*t744+t739.*t752];
mt2 = [t342+t354+t358+t361+t363+t371-t743.*(-t234+t256+t267+t283+t309+t330+t624+t650)+t658.*t671-t661.*t670+t676.*t684+t678.*t687+t694.*t700+t695.*t701+t710.*t716+t712.*t715+t727.*t734+t726.*t740-t739.*t746,t341+t347+t362+t374+t375+t376+t659.*t667-t660.*t666-t674.*t685-t675.*t686-t690.*t702-t691.*t703-t706.*t718-t708.*t719-t722.*t738-t723.*t742-t733.*t748-t735.*t749,et7+et8+et9+et10+et11+et12,t345+t349+t380+t383+t385+t388-t735.*(-t240+t258+t266+t308+t312+t328+t640+t653)-t659.*t673+t660.*t672+t674.*t688+t675.*t689+t690.*t704+t691.*t705+t706.*t721+t708.*t720-t722.*t745+t723.*t744-t733.*t752];
mt3 = [t343+t350+t356+t359+t365+t369+t735.*(-t234+t256+t267+t283+t309+t330+t624+t650)-t659.*t670-t660.*t671-t674.*t684-t675.*t687-t690.*t700-t691.*t701-t706.*t716-t708.*t715-t723.*t734-t722.*t740+t733.*t746,t344+t367+t378+t382+t384+t389+t749.*(-t129+t148+t155+t196+t201+t216+t506+t516)-t664.*t666-t665.*t667-t680.*t685-t681.*t686-t696.*t702-t697.*t703-t711.*t719-t713.*t718+t729.*t738-t728.*t742+t737.*t748];
mt4 = [t345+t349+t380+t383+t385+t388-t751.*(-t129+t148+t155+t196+t201+t216+t506+t516)+t664.*t669-t665.*t668+t680.*t682+t681.*t683+t696.*t698+t697.*t699+t713.*t714+t711.*t717-t729.*t732+t728.*t736-t737.*t750,et13+et14+et15+et16+et17+et18,t336+t351+t370+t372+t373+t379-(-t129+t148+t155+t196+t201+t216+t506+t516).*(-t234+t256+t267+t283+t309+t330+t624+t650)-t664.*t671+t665.*t670-t680.*t684-t681.*t687-t696.*t700-t697.*t701-t711.*t715-t713.*t716-t728.*t734+t729.*t740-t737.*t746];
mt5 = [t342+t354+t358+t361+t363+t371-t749.*(-t125+t146+t163+t181+t198+t218+t486+t512)-t662.*t667+t663.*t666+t677.*t685+t679.*t686+t692.*t702+t693.*t703+t707.*t719+t709.*t718+t725.*t738+t724.*t742-t730.*t748,t343+t350+t356+t359+t365+t369+t751.*(-t125+t146+t163+t181+t198+t218+t486+t512)-t662.*t668-t663.*t669-t677.*t682-t679.*t683-t692.*t698-t693.*t699-t709.*t714-t707.*t717-t725.*t732-t724.*t736+t730.*t750];
mt6 = [t336+t351+t370+t372+t373+t379-(-t125+t146+t163+t181+t198+t218+t486+t512).*(-t240+t258+t266+t308+t312+t328+t640+t653)+t662.*t673-t663.*t672-t677.*t688-t679.*t689-t692.*t704-t693.*t705-t707.*t720-t709.*t721-t724.*t744+t725.*t745-t730.*t752,et19+et20+et21+et22+et23+et24];
MMred = reshape([mt1,mt2,mt3,mt4,mt5,mt6],4,4);
end