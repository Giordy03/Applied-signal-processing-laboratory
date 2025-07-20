clc, clear, close all
% 4/4 --> BPM = 60
T = 1;
fs = 44100;

A5 = 880;
A4 = 440;
A3 = 220;
A2 = 110;
A1 = 55;

E5f = A5 * 2^(-6/12);
E5 = A5 * 2^(-5/12);
C5 = A5 * 2^(-9/12);
C5s = A5 * 2^(-8/12);
C5f = A5 * 2^(-10/12);
F5 = A5 * 2^(-4/12);
G5 = A5 * 2^(-2/12);
D5 = A5 * 2^(-7/12);
D5s = A5 * 2^(-6/12);
D5f = A5 * 2^(-8/12);
B5f = A5 * 2^(1/12);


E4f = A4 * 2^(-6/12);
E4 = A4 * 2^(-5/12);
C4f = A4 * 2^(-10/12);
C4 = A4 * 2^(-9/12);
C4s = A4 * 2^(-8/12);
F4 = A4 * 2^(-4/12);
F4f = A4 * 2^(-5/12);
F4s = A4 * 2^(-3/12);
G4 = A4 * 2^(-2/12);
G4f = A4 * 2^(-3/12);
D4 = A4 * 2^(-7/12);
D4s = A4 * 2^(-6/12);
D4f = A4 * 2^(-8/12);
B4f = A4 * 2^(1/12);
A4f = A4 * 2^(-1/12);

E3f = A3 * 2^(-6/12);
E3 = A3 * 2^(-5/12);
C3 = A3 * 2^(-9/12);
C3s = A2 * 2^(-8/12);
F3 = A3 * 2^(-4/12);
G3 = A3 * 2^(-2/12);
D3 = A3 * 2^(-7/12);
D3s = A3 * 2^(-6/12);
A3f = A3 * 2^(-1/12);
B3f = A3 * 2^(1/12);

E2f = A2 * 2^(-6/12);
E2 = A2 * 2^(-5/12);
C2 = A2 * 2^(-9/12);
C2s = A2 * 2^(-8/12);
F2 = A2 * 2^(-4/12);
G2 = A2 * 2^(-2/12);
D2 = A2 * 2^(-7/12);
D2s = A2 * 2^(-6/12);
B2f = A2 * 2^(1/12);
B2 = A2 * 2^(2/12);
A2f = A2 * 2^(-1/12);

E1f = A1 * 2^(-6/12);
E1 = A1 * 2^(-5/12);
C1 = A1 * 2^(-9/12);
F1 = A1 * 2^(-4/12);
G1 = A1 * 2^(-2/12);
D1 = A1 * 2^(-7/12);
D1s = A1 * 2^(-6/12);
B1f = A1 * 2^(1/12);
B1 = A1 * 2^(2/12);
A1f = A1 * 2^(-1/12);

d4 = add_notes(D4, T/2);
d4_3 = add_notes(D4, 3*T/2);
c4 = add_notes(C4, T/2);
c4_1 = add_notes(C4, T);
e3f = add_notes(E3f, T/2);
b3f = add_notes(B3f, T/2);
b3f_2 = add_notes(B3f, 2*T);
f3 = add_notes(F3, T/2);
a3 = add_notes(A3, T/2);
f3_3 = add_notes(F3, 3*T/2);
g = add_notes([G1 G2], 4*T);
bf = add_notes([B1f B2f], 4*T);
ef = add_notes([E2f E3f], 4*T);
c = add_notes([C2 C3], 4*T);
f = add_notes([F2 F3], 4*T);
c4f_1 = add_notes(C4s, T);
b3f_1 = add_notes(B3f, T);
a3_1 = add_notes(A3, T);
ef_2 = add_notes([E2f E3f], 2*T);
d_2 = add_notes([D2 D3], 2*T);
cs_2 = add_notes([C2s C3s], 2*T);
bf_2 = add_notes([B1f B2f], 2*T);
g_1 = add_notes([G3 G4], T);
f_1 = add_notes([F3 F4], T);
c_2 = add_notes([C2 C3], 2*T);
a_2 = add_notes([A1 A2], 2*T);
b3f_s = add_notes(B3f, T/4);
b3f_l = add_notes(B3f, 7*T/4);
b_1 = add_notes([B1 B2], T);
bf_1 = add_notes([B1f B2f], T);
af_1 = add_notes([A1f A2f], T);
glow_1 = add_notes([G1 G2], T);
ef_3t = add_notes([E3f E2f], 3*T);
d_1 = add_notes([D3 D2], T);
f_3 = add_notes([F2 F3], 3/2*T);
e = add_notes([E2 E3], T/2);
ef_1 = add_notes([E2f E4f], T);
bf_3t = add_notes([B1f B2f], 3*T);
af_2 = add_notes([A1f A2f], 2*T);
ef_6 = add_notes([E2f E3f], 6*T);
bf_l = add_notes([B1f B2f], 7/2*T);


left_1 = [d4 d4 d4 d4 d4 d4_3 c4 c4 d4 c4 b3f_s b3f_l c4 c4 c4 c4_1 c4 f3 f3 b3f b3f c4 b3f a3 f3_3];
left_2 = [g bf ef c];
left_3 = [f c4f_1 b3f_1 a3_1 b3f_1 c4f_1 b3f_1 a3_1 b3f_1 ef_2 d_2]; 
left_4 = [cs_2 cs_2 c bf_2 g_1 f_1 bf_2 g_1 f_1];
left_5 = [f g c c_2 a_2];
left_6 = [bf g c_2 b_1 bf_1 a_2 af_1 glow_1];
left_7 = [ef_3t d_1 c f_3 e ef_1 d_1 bf_3t bf_1];
left_8 = [ef_2 d_2 c_2 af_2 ef_6];
left_9 = [bf_2 g_1 f_1 bf_2 g_1 f_1 bf_l];

%% RIGHT
n1 = add_notes([F4 G4 B4f], T/2);
n1_3 = add_notes([F4 G4 B4f], 3*T/2);
n2 = add_notes([E4 G4 B4f], T/2);
n2f = add_notes([E4f G4 B4f], T/2);
n3_s = add_notes([D4 F4 G4], T/4);
n4_l = add_notes([C4 E4 G4], 7*T/4);
n5 = add_notes([E4f F4 A4], T/2);
n6_1 = add_notes([E4f G4 B4f], T);
f4 = add_notes(F4, T/2);
n7 = add_notes([D4 F4 B4f], T/2);
n8_3 = add_notes([B3f D4 F4], 3*T/2);
n_new = add_notes([D4 F4 A4], T/2);
n9 = add_notes([B3f D4], T/2);
n9_2 = add_notes([B3f D4], 2*T);
n10 = add_notes([A3f D4], T/2);
n11 = add_notes([A3f D4 E4f], T/2);
n11_1 = add_notes([A3f D4 E4f], T);
b3f_1 = add_notes(B3f, T);
g4 = add_notes(G4, T/2);
g4_1 = add_notes(G4, T);
b4f = add_notes(B4f, T/2);
e4f = add_notes(E4f, T/2);
n12 = add_notes([A4f D5], T/2);
e5f = add_notes(E5f, T/2);
n13 = add_notes([G4 B4f], T/2);
n14 = add_notes([E4f, G4], T/2);
n14_1 = add_notes([E4f, G4], T);
n14_3 = add_notes([E4f, G4], 3*T/2);
n15 = add_notes([A3, F4], T/2);
n16 = add_notes([B3f, G4], T/2);
n17 = add_notes([C4, A4], T/2);
n18 = add_notes([D4, B4f], T/2);
n19 = add_notes([E4f G4f C5f], T/2);
n20 = add_notes([D4 F4 B4f], T/2);
n20_1 = add_notes([D4 F4 B4f], T);
n21 = add_notes([C4 E4 A4], T/2);
n21 = add_notes([C4s E4 A4], T/2);
g4_s = add_notes(G4, T/4);
f4_l  = add_notes(F4, 5*T/4);
b4f_1 = add_notes(B4f, T);
e4 = add_notes(E4, T/2);
f4_1 = add_notes(F4, T);
a4_3t = add_notes(A4, 3*T);
a4_3 = add_notes(A4, 3*T/2);
d4_3t = add_notes(D4, 3*T);
d4_s = add_notes(D4, T/4);
d4_ll = add_notes(D4, 13*T/4);
d4_l = add_notes(D4, 9*T/4);
pause = zeros(1, T/2*fs);
f4_s = add_notes(F4, T/4);
c4_l = add_notes(C4, 5*T/4);
e4f_s = add_notes(E4f, T/4);
c4_ll = add_notes(C4, 9*T/4);
% b3f = add_notes(B3f, T/2);
c4_s = add_notes(C4, T/4);
e4f_s = add_notes(E4f, T/4);
d4_3 = add_notes(D4, 3*T/2);
pause_l = zeros(1, 3/4*T*fs);
g4_l = add_notes(G4, 3*T/4);
a4_1 = add_notes(A4, T);
g4_2 = add_notes(G4, 2*T);
e4f_1 = add_notes(E4f, T);
d4_1 = add_notes(D4, T);
g4_3t = add_notes(G4, 3*T);
a4f_s = add_notes(A4f, T/4);
a4f_3 = add_notes(A4f, 3*T/2);
a4f = add_notes(A4f, T/2);
b4f_s = add_notes(B4f, T/4);
d4_lm = add_notes(D4, 11*T/4);
a4_m = add_notes(A4, 3*T/4);
b4f_m = add_notes(B4f, 3*T/4);
c4_4 = add_notes(C4, 4*T);
g4_lm = add_notes(G4, 11*T/4);
a4f_m = add_notes(A4f, 3*T/4);
d4_m = add_notes(D4, 5/2*T);
a4f_s = add_notes(A4f, T/4);
e4f_3 = add_notes(E4f, 3*T/2);
c4f = add_notes(C4f, T/2);
d4f = add_notes(D4f, T/2);
c4f_s = add_notes(C4f, T/4);
d4f_s = add_notes(D4f, T/4);
s3f_m = add_notes(B3f, 3/4*T);
n22 = add_notes([A4f, C5], T/2);
n23 = add_notes([F4s A4], T/2);
n24 = add_notes([F4 A4f], T/2);
a4 = add_notes(A4, T/2);

right_1 = [n1 n1 n1 n1 n1 n1_3 n2 n2 n1 n2 n3_s n4_l n5 n5 n5 n6_1 n5 f4 f4 n7 n7 n2f n7 n_new n8_3];
right_2 = [n9 n9 n9 n9_2 n9 n10 n10 n11 n11_1 b3f_1 g4_1 e4f g4 b4f n12 e5f n13 e5f n14 n14 n14 n14_1 n14_3];
right_3 = [n15 n15 g4 n15 e4f_s c4_s n16 n17 n18 n19 n19 n20_1 n21 n21 n20_1 n19 n19 n20_1 n21 n21 n20_1 g4 g4 g4 g4_s f4_l b4f_1];
right_4 = [e4 e4 e4 e4 f4_1 f3 a3 a4_3t a3 b3f b3f f3 b3f d4_m b3f f3 b3f d4_m];
right_5 = [d4_s d4_ll b3f c4 d4_s d4_l pause d4_s d4_s e4f f4_s e4f d4 c4_l c4 d4 e4f_s f4 e4f d4 c4_ll];
right_6 = [d4_s d4_lm d4 f4 a4_m g4_s g4_2 pause g4 b4f b4f b4f b4f b4f_m g4_s e4f_1 c4_4];
right_7 = [g4_s g4_lm f4 g4_s a4f_s g4_3t pause g4 a4f_m g4_s g4 f4 f4_1 pause_l b3f_s b3f f4 f4 g4_s g4_l a4f a4f b4f_s a4f_s];
right_8 = [a4 g4_1 f4_s g4_s e4f_3 f4_s g4_s e4f_3 b3f_s b3f_s c4f d4f c4f_s d4f_s c4f_s s3f_m e4f g4 b4f n22 e5f n22 e5f n23 e5f n24 e5f];
right_9 = [b3f f3 b3f d4_m b3f f3 b3f d4_m d4_s d4_ll];
music = [right_1+left_1, right_2+left_2, right_3+left_3, right_4+left_4, right_5+left_5, right_6+left_6, right_7+left_7, right_8+left_8, right_9+left_9];
music_obj = audioplayer(music, fs);
playblocking(music_obj)
audiowrite("bohemian_rhapsody.wav", music, fs);

function smooth_sound = add_notes(freqs, T)
    fs = 44100;
    t = 0:1/fs:T-1/fs;
    notes = zeros(T*fs, length(freqs));
    for i = 1:length(freqs)
        notes(:, i) = note_generation(freqs(i), T);
    end
    note = notes/max(abs(notes));
    ADSR = interp1([0 0.01 0.9 1]*T, [0 1 0.1 0], t, 'pchip');
    smooth_sound = note.'.*ADSR;
end

function note = note_generation(f0, T)
    fs = 44100;
    t = 0:1/fs:T-1/fs;
    init_tone = (sin(2*pi*f0*t) ...
        + sin(2*pi*f0*t + 0.2/f0) + sin(2*pi*f0*t - 0.2/f0) ...
        + sin(2*pi*f0*t + 0.4/f0) + sin(2*pi*f0*t - 0.4/f0) ...
        + sawtooth(2*pi*f0*t + 0.3/f0) + sawtooth(2*pi*f0*t - 0.3/f0) ...
        + sawtooth(2*pi*f0*t + 0.5/f0) + sawtooth(2*pi*f0*t - 0.5/f0) ...
        );
   
    lfo = f0/16;
    lfo_sine = sin(2*pi*lfo*t).*(-t/(20*T));

    lfo_2 = f0/8;
    lfo_sine_2 = 0.12*sin(2*pi*lfo_2*t)+ 1;
    
    pwm = init_tone < lfo_sine;
    pwm = pwm .* lfo_sine_2;

    Bt_lp = 7*f0*2*pi/fs;
    Wp_lp = f0*2*pi/fs;
    Ws_lp = 8*f0*2*pi/fs;
    
    Wc_lp = (Wp_lp+Ws_lp)/2;

    N_lp = ceil((6.1*pi)/Bt_lp);
    n_lp = 0:N_lp-1;

    window_lp = bartlett(N_lp);
    window_lp = window_lp';

    M_lp = ceil((N_lp-1)/2);

    h_id_lp=Wc_lp/pi*sinc(Wc_lp/pi*(n_lp-M_lp));
    
    h_lp = window_lp .* h_id_lp;
    h_lp = h_lp/sum(h_lp);

    Bt_hp = 100*2*pi/fs;
    Wp_hp = 60*2*pi/fs;
    Ws_hp = 20*2*pi/fs;
    Wc_hp = (Wp_hp+Ws_hp)/2;

    N_hp = ceil((6.1 * pi) / Bt_hp);
    n_hp = 0 : N_hp - 1;

    window_hp = bartlett(N_hp);
    window_hp = window_hp';

    M_hp = ceil((N_hp-1)/2);
    
    delta_hp = zeros(1, N_hp);
    delta_hp(M_hp+1) = 1;
    h_id_hp = delta_hp - (Wc_hp/pi * sinc(Wc_hp/pi * (n_hp - M_hp)));
    
    h_hp = window_hp .* h_id_hp;
    h_hp = h_hp/sum(h_hp);
    
    note = filter(h_lp, 1, pwm);
    note = filter(h_hp, 1, note);
end
