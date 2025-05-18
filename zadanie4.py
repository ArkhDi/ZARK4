import numpy as np
from scipy.special import spherical_jn, spherical_yn
import requests
import matplotlib.pyplot as plt

class ESACalculator:
    def __init__(self, diameter, fmin, fmax):
        self.diameter = diameter
        self.radius = diameter / 2
        self.fmin = fmin
        self.fmax = fmax

    def calculate_esa(self, frequency):
        k = 2 * np.pi * frequency / 3e8  
        r = self.radius
        lambda_ = 3e8 / frequency  

        n_max = 10
        rcs = (lambda_**2 / np.pi) * abs(sum(((-1)**n) * (n + 0.5) * (self.b_n(k, r, n) - self.a_n(k, r, n)) for n in range(1, n_max + 1)))**2

        return rcs

    def a_n(self, k, r, n):
        jn = spherical_jn(n , k * r)
        hn = spherical_jn(n, k * r) + 1j * spherical_yn(n, k * r)
        return jn / hn

    def b_n(self, k, r, n):
        jn_minus_1=spherical_jn(n-1, k * r)
        hn_minus_1=spherical_jn(n , k * r) + 1j * spherical_yn(n, k * r)
        jn = spherical_jn(n , k * r)
        hn = spherical_jn(n , k * r) + 1j * spherical_yn(n, k * r)
        
        
        return (k * r * jn_minus_1 - n * jn) / (k * r * hn_minus_1 - n * hn)

class ESAResultSaver:
    def __init__(self, filename):
        self.filename = filename

    def save_results(self, frequencies, esa_values):
        with open(self.filename, 'w') as file:
            for freq, rcs in zip(frequencies, esa_values):
                file.write(f"{freq:14.6e}{'':4}{esa:14.6e}\n")
   

def main():
    url = "https://jenyay.net/uploads/Student/Modelling/task_rcs_01.txt"
    response = requests.get(url)
    lines = response.text.splitlines()
    diameter_str = lines[1].split('=')[1].split(';')[0].strip()
    fmin_str = lines[1].split('=')[2].split(';')[0].strip()
    fmax_str = lines[1].split('=')[3].split(';')[0].strip()

    diameter = float(diameter_str)
    print(diameter)
    fmin = float(fmin_str)
    print(fmin)
    fmax = float(fmax_str)
    print(fmax)

   
    esa_calculator = ESACalculator(diameter, fmin, fmax)

    frequencies = np.linspace(fmin, fmax, 400)

    esa_values = [esa_calculator.calculate_esa(freq) for freq in frequencies]

    
    result_saver = ESAResultSaver("results.txt")
    result_saver.save_results(frequencies, esa_values)
    
    plt.figure(figsize=(10, 10)) 
    plt.plot(frequencies, esa_values)
    plt.xlabel("FREQ (Hz)")
    plt.ylabel("ESA (m²)")
    plt.title("ESA FROM FREQ")
    plt.grid(True)  
    plt.show()
   

if __name__ == "__main__":
    main()
