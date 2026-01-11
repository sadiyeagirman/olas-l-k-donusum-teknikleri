import sympy as sp

# Sembolik işlemler için değişkenler
x, y, t = sp.symbols('x y t', real=True)

# 1. d.f. (Dağılım Fonksiyonu) Tekniği
# Bu fonksiyon F_Y(y) = P(g(X) <= y) eşitliğini kulanarak çözüm yapar.
def teknik_df(g_x, f_x=None, F_x=None, alt=0, ust=1):
    try:
        # F_x yoksa f_x üzerinden integral ile hesaplanır
        if F_x is None and f_x is not None:
            Fx = sp.integrate(f_x, (x, alt, x))
        else:
            Fx = F_x
            
        # g(x)=y denkleminden x çekilir (ters fonksiyon)
        ters_x = sp.solve(sp.Eq(y, g_x), x)[0]
        
        # Fonksiyonun artan veya azalan olma durumu kontrol edilir
        egim = sp.diff(g_x, x).subs(x, (alt + ust) / 2)
        
        if egim > 0:
            Fy = Fx.subs(x, ters_x) # g artan ise
        else:
            Fy = 1 - Fx.subs(x, ters_x) # g azalan ise
            
        return sp.simplify(Fy)
    except:
        return "d.f. tekniği ile sonuç bulunamadı."

# 2. o.y.f. (Olasılık Yoğunluk Fonksiyonu) Tekniği
# Bu fonksiyon f_Y(y) = f_X(g^-1(y)) * |d(g^-1(y))/dy| formülünü uygular.
def teknik_oyf(g_x, f_x, alt=0, ust=1):
    try:
        # g(x)=y denkleminin kökleri (ters fonksiyonlar) bulunur
        kokler = sp.solve(sp.Eq(y, g_x), x)
        toplam_oyf = 0
        
        for k in kokler:
            # Her kök için mutlak değerli türev çarpan hesaplanır
            turev_carpan = sp.abs(sp.diff(k, y))
            toplam_oyf += f_x.subs(x, k) * turev_carpan
            
        return sp.simplify(toplam_oyf)
    except:
        return "o.y.f. tekniği ile sonuç bulunamadı."

# 3. m.ç.f. (Moment Çıkaran Fonksiyon) Tekniği
# Bu fonksiyon M_Y(t) = E[e^(t*g(X))] beklenen değerini hesaplar.
def teknik_mcf(g_x, f_x, alt=0, ust=1):
    try:
        # Dönüşümün moment fonksiyonu integrali alınır
        m_y = sp.integrate(sp.exp(t * g_x) * f_x, (x, alt, ust))
        return sp.simplify(m_y)
    except:
        return "m.ç.f. tekniği ile sonuç bulunamdı."
