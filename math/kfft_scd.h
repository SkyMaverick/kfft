#define K1 48

/* Frac: table of 2^(-k), k=1,...,K1 */
static double Frac[K1] = {
    5.0000000000000000e-001, /* 2^{-1} */
    2.5000000000000000e-001, /* 2^{-2} */
    1.2500000000000000e-001, /* 2^{-3} */
    6.2500000000000000e-002, /* 2^{-4} */
    3.1250000000000000e-002, /* 2^{-5} */
    1.5625000000000000e-002, /* 2^{-6} */
    7.8125000000000000e-003, /* 2^{-7} */
    3.9062500000000000e-003, /* 2^{-8} */
    1.9531250000000000e-003, /* 2^{-9} */
    9.7656250000000000e-004, /* 2^{-10} */
    4.8828125000000000e-004, /* 2^{-11} */
    2.4414062500000000e-004, /* 2^{-12} */
    1.2207031250000000e-004, /* 2^{-13} */
    6.1035156250000000e-005, /* 2^{-14} */
    3.0517578125000000e-005, /* 2^{-15} */
    1.5258789062500000e-005, /* 2^{-16} */
    7.6293945312500000e-006, /* 2^{-17} */
    3.8146972656250000e-006, /* 2^{-18} */
    1.9073486328125000e-006, /* 2^{-19} */
    9.5367431640625000e-007, /* 2^{-20} */
    4.7683715820312500e-007, /* 2^{-21} */
    2.3841857910156250e-007, /* 2^{-22} */
    1.1920928955078130e-007, /* 2^{-23} */
    5.9604644775390630e-008, /* 2^{-24} */
    2.9802322387695310e-008, /* 2^{-25} */
    1.4901161193847660e-008, /* 2^{-26} */
    7.4505805969238280e-009, /* 2^{-27} */
    3.7252902984619140e-009, /* 2^{-28} */
    1.8626451492309570e-009, /* 2^{-29} */
    9.3132257461547850e-010, /* 2^{-30} */
    4.6566128730773930e-010, /* 2^{-31} */
    2.3283064365386960e-010, /* 2^{-32} */
    1.1641532182693480e-010, /* 2^{-33} */
    5.8207660913467410e-011, /* 2^{-34} */
    2.9103830456733700e-011, /* 2^{-35} */
    1.4551915228366850e-011, /* 2^{-36} */
    7.2759576141834260e-012, /* 2^{-37} */
    3.6379788070917130e-012, /* 2^{-38} */
    1.8189894035458560e-012, /* 2^{-39} */
    9.0949470177292820e-013, /* 2^{-40} */
    4.5474735088646410e-013, /* 2^{-41} */
    2.2737367544323210e-013, /* 2^{-42} */
    1.1368683772161600e-013, /* 2^{-43} */
    5.6843418860808010e-014, /* 2^{-44} */
    2.8421709430404010e-014, /* 2^{-45} */
    1.4210854715202000e-014, /* 2^{-46} */
    7.1054273576010020e-015, /* 2^{-47} */
    3.5527136788005010e-015, /* 2^{-48} */
};

/* Cosi: table of cos(2*pi*2^(-k)) k=1,...,K1 */
static double Cosi[K1] = {
    -1.0000000000000000e+000, /* cos(2*PI*2^{-1}) */
    0.0000000000000000e+000,  /* cos(2*PI*2^{-2}) */
    7.0710678118654760e-001,  /* cos(2*PI*2^{-3}) */
    9.2387953251128670e-001,  /* cos(2*PI*2^{-4}) */
    9.8078528040323040e-001,  /* cos(2*PI*2^{-5}) */
    9.9518472667219690e-001,  /* cos(2*PI*2^{-6}) */
    9.9879545620517240e-001,  /* cos(2*PI*2^{-7}) */
    9.9969881869620420e-001,  /* cos(2*PI*2^{-8}) */
    9.9992470183914450e-001,  /* cos(2*PI*2^{-9}) */
    9.9998117528260110e-001,  /* cos(2*PI*2^{-10}) */
    9.9999529380957620e-001,  /* cos(2*PI*2^{-11}) */
    9.9999882345170190e-001,  /* cos(2*PI*2^{-12}) */
    9.9999970586288220e-001,  /* cos(2*PI*2^{-13}) */
    9.9999992646571790e-001,  /* cos(2*PI*2^{-14}) */
    9.9999998161642930e-001,  /* cos(2*PI*2^{-15}) */
    9.9999999540410730e-001,  /* cos(2*PI*2^{-16}) */
    9.9999999885102690e-001,  /* cos(2*PI*2^{-17}) */
    9.9999999971275670e-001,  /* cos(2*PI*2^{-18}) */
    9.9999999992818920e-001,  /* cos(2*PI*2^{-19}) */
    9.9999999998204720e-001,  /* cos(2*PI*2^{-20}) */
    9.9999999999551180e-001,  /* cos(2*PI*2^{-21}) */
    9.9999999999887800e-001,  /* cos(2*PI*2^{-22}) */
    9.9999999999971940e-001,  /* cos(2*PI*2^{-23}) */
    9.9999999999992980e-001,  /* cos(2*PI*2^{-24}) */
    9.9999999999998250e-001,  /* cos(2*PI*2^{-25}) */
    9.9999999999999570e-001,  /* cos(2*PI*2^{-26}) */
    9.9999999999999890e-001,  /* cos(2*PI*2^{-27}) */
    9.9999999999999980e-001,  /* cos(2*PI*2^{-28}) */
    9.9999999999999990e-001,  /* cos(2*PI*2^{-29}) */
    1.0000000000000000e+000,  /* cos(2*PI*2^{-30}) */
    1.0000000000000000e+000,  /* cos(2*PI*2^{-31}) */
    1.0000000000000000e+000,  /* cos(2*PI*2^{-32}) */
    1.0000000000000000e+000,  /* cos(2*PI*2^{-33}) */
    1.0000000000000000e+000,  /* cos(2*PI*2^{-34}) */
    1.0000000000000000e+000,  /* cos(2*PI*2^{-35}) */
    1.0000000000000000e+000,  /* cos(2*PI*2^{-36}) */
    1.0000000000000000e+000,  /* cos(2*PI*2^{-37}) */
    1.0000000000000000e+000,  /* cos(2*PI*2^{-38}) */
    1.0000000000000000e+000,  /* cos(2*PI*2^{-39}) */
    1.0000000000000000e+000,  /* cos(2*PI*2^{-40}) */
    1.0000000000000000e+000,  /* cos(2*PI*2^{-41}) */
    1.0000000000000000e+000,  /* cos(2*PI*2^{-42}) */
    1.0000000000000000e+000,  /* cos(2*PI*2^{-43}) */
    1.0000000000000000e+000,  /* cos(2*PI*2^{-44}) */
    1.0000000000000000e+000,  /* cos(2*PI*2^{-45}) */
    1.0000000000000000e+000,  /* cos(2*PI*2^{-46}) */
    1.0000000000000000e+000,  /* cos(2*PI*2^{-47}) */
    1.0000000000000000e+000,  /* cos(2*PI*2^{-48}) */
};

/* Sine: table of sin(2*pi*2^(-k)) k=1,...,K1 */
static double Sine[K1] = {
    0.0000000000000000e-000, /* sin(2*Pi*2^{-1}) */
    1.0000000000000000e+000, /* sin(2*Pi*2^{-2}) */
    7.0710678118654750e-001, /* sin(2*Pi*2^{-3}) */
    3.8268343236508980e-001, /* sin(2*Pi*2^{-4}) */
    1.9509032201612820e-001, /* sin(2*Pi*2^{-5}) */
    9.8017140329560600e-002, /* sin(2*Pi*2^{-6}) */
    4.9067674327418010e-002, /* sin(2*Pi*2^{-7}) */
    2.4541228522912290e-002, /* sin(2*Pi*2^{-8}) */
    1.2271538285719930e-002, /* sin(2*Pi*2^{-9}) */
    6.1358846491544750e-003, /* sin(2*Pi*2^{-10}) */
    3.0679567629659760e-003, /* sin(2*Pi*2^{-11}) */
    1.5339801862847660e-003, /* sin(2*Pi*2^{-12}) */
    7.6699031874270450e-004, /* sin(2*Pi*2^{-13}) */
    3.8349518757139560e-004, /* sin(2*Pi*2^{-14}) */
    1.9174759731070330e-004, /* sin(2*Pi*2^{-15}) */
    9.5873799095977340e-005, /* sin(2*Pi*2^{-16}) */
    4.7936899603066880e-005, /* sin(2*Pi*2^{-17}) */
    2.3968449808418220e-005, /* sin(2*Pi*2^{-18}) */
    1.1984224905069710e-005, /* sin(2*Pi*2^{-19}) */
    5.9921124526424280e-006, /* sin(2*Pi*2^{-20}) */
    2.9960562263346610e-006, /* sin(2*Pi*2^{-21}) */
    1.4980281131690110e-006, /* sin(2*Pi*2^{-22}) */
    7.4901405658471570e-007, /* sin(2*Pi*2^{-23}) */
    3.7450702829238410e-007, /* sin(2*Pi*2^{-24}) */
    1.8725351414619530e-007, /* sin(2*Pi*2^{-25}) */
    9.3626757073098080e-008, /* sin(2*Pi*2^{-26}) */
    4.6813378536549090e-008, /* sin(2*Pi*2^{-27}) */
    2.3406689268274550e-008, /* sin(2*Pi*2^{-28}) */
    1.1703344634137280e-008, /* sin(2*Pi*2^{-29}) */
    5.8516723170686380e-009, /* sin(2*Pi*2^{-30}) */
    2.9258361585343190e-009, /* sin(2*Pi*2^{-31}) */
    1.4629180792671600e-009, /* sin(2*Pi*2^{-32}) */
    7.3145903963357980e-010, /* sin(2*Pi*2^{-33}) */
    3.6572951981678990e-010, /* sin(2*Pi*2^{-34}) */
    1.8286475990839500e-010, /* sin(2*Pi*2^{-35}) */
    9.1432379954197480e-011, /* sin(2*Pi*2^{-36}) */
    4.5716189977098740e-011, /* sin(2*Pi*2^{-37}) */
    2.2858094988549370e-011, /* sin(2*Pi*2^{-38}) */
    1.1429047494274680e-011, /* sin(2*Pi*2^{-39}) */
    5.7145237471373420e-012, /* sin(2*Pi*2^{-40}) */
    2.8572618735686710e-012, /* sin(2*Pi*2^{-41}) */
    1.4286309367843360e-012, /* sin(2*Pi*2^{-42}) */
    7.1431546839216780e-013, /* sin(2*Pi*2^{-43}) */
    3.5715773419608390e-013, /* sin(2*Pi*2^{-44}) */
    1.7857886709804190e-013, /* sin(2*Pi*2^{-45}) */
    8.9289433549020970e-014, /* sin(2*Pi*2^{-46}) */
    4.4644716774510490e-014, /* sin(2*Pi*2^{-47}) */
    2.2322358387255240e-014, /* sin(2*Pi*2^{-48}) */
};
