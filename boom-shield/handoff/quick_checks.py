"""Quick back-of-envelope checks behind the 2026-10-20 review addendum.

Standalone (numpy only). Prints:
  1. Air and helium properties at 20 km (US Standard Atmosphere 1976) and the
     ideal buoyant capacity of a 15 m elliptical helium envelope.
  2. Limp-panel areal mass needed for a given pressure-amplitude transmission A
     at normal incidence (mass law, air impedance at 20 km).
  3. Single limp sheet vs double wall (two membranes + 15 m helium gap)
     transmission versus frequency, and the peak transmitted fraction of an
     N-wave of duration tau.
  4. Helium inventory for a 4,000 km corridor vs world production (~170e6 m3/yr STP).
  5. Rough boom overpressure 250 m below a Mach 3.2 aircraft (Carlson, Kr = 1)
     compared with deck dead load.
  6. Tether drag through the jet stream.
  7. Mass-law vs buoyancy closure: peak transmission when ALL lift is spent on
     acoustically active mass, for envelope thickness T.

These are screening estimates, not designs. See boom-shield/sim/ for the
propagation model that supersedes them.
"""
import numpy as np

# --- Atmosphere at 20 km (US Std Atm 1976) ---------------------------------
T = 216.65          # K
P = 5474.9          # Pa
R_AIR = 287.05
R_HE = 2077.1
G = 9.80665

rho = P / (R_AIR * T)
c = np.sqrt(1.4 * R_AIR * T)
Z = rho * c
rho_he = P / (R_HE * T)
c_he = np.sqrt(5 / 3 * R_HE * T)
Z_he = rho_he * c_he


def print_atmosphere():
    print(f"air rho={rho:.4f} c={c:.1f} Z={Z:.1f}; He rho={rho_he:.4f} c={c_he:.0f} Z={Z_he:.1f}")
    print("buoyancy cap T=15:", (rho - rho_he) * np.pi * 15 / 4, " vacuum:", rho * np.pi * 15 / 4)
    print("T for 4 kg/m2:", 4 * 4 / (np.pi * (rho - rho_he)))


def print_limp_mass_needed():
    for f in (5, 10, 20):
        for A in (0.4, 0.2, 0.08, 0.05):
            r = np.sqrt(1 / A**2 - 1)
            m = r * 2 * Z / (2 * np.pi * f)
            print(f"f={f}Hz A={A}: limp mass {m:.1f} kg/m2")


# --- Transfer functions (normal incidence) ----------------------------------
def H_single(m):
    """Limp sheet of areal mass m between air half-spaces."""
    return lambda f: 1 / (1 + 1j * 2 * np.pi * f * m / (2 * Z))


def H_double(m1, m2, d):
    """Membrane m1 | helium layer d | membrane m2, air on both sides."""
    def h(f):
        w = 2 * np.pi * f
        k = w / c_he
        M = lambda m: np.array([[1, 1j * w * m], [0, 1]])
        L = np.array([[np.cos(k * d), 1j * Z_he * np.sin(k * d)],
                      [1j * np.sin(k * d) / Z_he, np.cos(k * d)]])
        A = M(m1) @ L @ M(m2)
        return 2 / (A[0, 0] + A[0, 1] / Z + A[1, 0] * Z + A[1, 1])
    return h


def peak_through(tau, H, fs=20000, Tw=4.0):
    """Peak |p_out| / peak p_in for a unit N-wave of duration tau."""
    n = int(Tw * fs)
    t = np.arange(n) / fs
    p = np.where((t > 0.5) & (t < 0.5 + tau), 1 - 2 * (t - 0.5) / tau, 0.0)
    F = np.fft.rfft(p)
    f = np.fft.rfftfreq(n, 1 / fs)
    Hf = np.array([H(fi) if fi > 0 else 1.0 for fi in f])
    return np.max(np.abs(np.fft.irfft(F * Hf, n)))


def print_transmission():
    print("\nAmplitude transmission |t|: single limp vs double wall (He gap 15 m)")
    for f in (2, 3, 5, 10, 20, 30, 50, 100):
        print(f"f={f:4}Hz single0.5={abs(H_single(0.5)(f)):.3f} "
              f"double0.25+0.25={abs(H_double(0.25, 0.25, 15)(f)):.3f} "
              f"double0.5+0.5={abs(H_double(0.5, 0.5, 15)(f)):.3f}")
    print("\nPeak transmitted / incident for an N-wave (linear, normal incidence, no propagation after)")
    for tau in (0.05, 0.1, 0.2):
        print(f"tau={tau}s: limp0.5={peak_through(tau, H_single(0.5)):.2f} "
              f"limp1={peak_through(tau, H_single(1)):.2f} "
              f"limp4={peak_through(tau, H_single(4)):.2f} "
              f"dbl0.25x2/15m={peak_through(tau, H_double(.25, .25, 15)):.2f} "
              f"dbl0.45x2/15m={peak_through(tau, H_double(.45, .45, 15)):.2f}")


def print_helium_inventory():
    for W, Tk in ((1000, 15), (250, 15), (250, 66)):
        V = np.pi * W * Tk / 4 * 4.0e6   # 4000 km corridor
        stp = V * P / 101325 * 273.15 / T
        print(f"W={W} T={Tk}: {V:.2e} m3 at altitude, {stp:.2e} m3 STP, "
              f"= {stp / 1.7e8:.1f} yr of ~170e6 m3/yr world He")


def print_boom_load():
    # Carlson simplified method, no ground reflection, h = 250 m below Mach 3.2, l = 32.7 m
    M, h, l = 3.2, 250, 32.7
    for Kp in (0.06, 0.1, 0.13):
        dp = Kp * P * (M**2 - 1)**0.125 * h**-0.75 * l**0.75
        print(f"Kp={Kp}: dp~{dp:.0f} Pa vs dead load 4kg/m2={4 * G:.0f} Pa, 0.9kg/m2={0.9 * G:.0f} Pa")


def print_tether_drag():
    for rho_t, U, d in ((0.41, 60, 0.05), (0.41, 80, 0.05)):
        q = 0.5 * rho_t * U**2 * 1.2 * d
        print("tether drag N/m:", q, " over 8 km:", q * 8000 / 1e3, "kN")


# --- Closure: all lift spent on acoustic mass -------------------------------
def _nwave_spectrum(tau, fs=8000, Tw=6.0):
    n = int(Tw * fs)
    t = np.arange(n) / fs
    p = np.where((t > 0.5) & (t < 0.5 + tau), 1 - 2 * (t - 0.5) / tau, 0.0)
    f = np.fft.rfftfreq(n, 1 / fs)
    f[0] = 1e-9
    return 2 * np.pi * f, np.fft.rfft(p), n


def closure_single(T_env, tau):
    w, F, n = _nwave_spectrum(tau)
    cap = (rho - rho_he) * np.pi * T_env / 4
    H = 1 / (1 + 1j * w * cap / (2 * Z))
    return np.max(np.abs(np.fft.irfft(F * H, n))), cap


def closure_double(T_env, tau):
    w, F, n = _nwave_spectrum(tau)
    cap = (rho - rho_he) * np.pi * T_env / 4
    m = cap / 2
    k = w / c_he
    d = T_env  # gas gap ~ envelope thickness
    a11 = np.cos(k * d) - w * m * np.sin(k * d) / Z_he
    a12 = 1j * w * m * np.cos(k * d) + 1j * Z_he * np.sin(k * d)
    a21 = 1j * np.sin(k * d) / Z_he
    a22 = np.cos(k * d)
    b11, b12 = a11, a11 * 1j * w * m + a12
    b21, b22 = a21, a21 * 1j * w * m + a22
    H = 2 / (b11 + b12 / Z + b21 * Z + b22)
    return np.max(np.abs(np.fft.irfft(F * H, n))), cap


def print_closure():
    print("All lift spent on acoustic mass (zero structure). N-wave tau=0.1 s, 20 km")
    for T_env in (15, 30, 60, 100, 200):
        s, cap = closure_single(T_env, 0.1)
        d, _ = closure_double(T_env, 0.1)
        print(f"T={T_env:4} m cap={cap:5.2f} kg/m2  single-limp peak A={s:.2f}  double-wall(He gap=T) A={d:.2f}")


if __name__ == "__main__":
    print_atmosphere()
    print_limp_mass_needed()
    print_transmission()
    print_helium_inventory()
    print_boom_load()
    print_tether_drag()
    print_closure()
