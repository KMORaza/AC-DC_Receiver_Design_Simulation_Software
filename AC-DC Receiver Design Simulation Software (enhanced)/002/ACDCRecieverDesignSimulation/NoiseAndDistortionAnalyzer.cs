using System;

namespace ACDCSimulator
{
    public class NoiseAndDistortionAnalyzer
    {
        private double inputVoltage;
        private double loadResistance;
        private string rectifierType;
        private double filterCapacitance;
        private double filterInductance;
        private double filterCutoff;
        private double pwmCarrierFreq;
        private double outputVoltage;

        public NoiseAndDistortionAnalyzer(
            double inputVoltage,
            double loadResistance,
            string rectifierType,
            double filterCapacitance,
            double filterInductance,
            double filterCutoff,
            double pwmCarrierFreq,
            double outputVoltage)
        {
            this.inputVoltage = inputVoltage;
            this.loadResistance = loadResistance;
            this.rectifierType = rectifierType;
            this.filterCapacitance = filterCapacitance; // in µF
            this.filterInductance = filterInductance;   // in mH
            this.filterCutoff = filterCutoff;           // in Hz
            this.pwmCarrierFreq = pwmCarrierFreq;       // in Hz
            this.outputVoltage = outputVoltage;         // in V
        }

        public double CalculateSNR()
        {
            // Signal power (based on output voltage)
            double signalPower = Math.Pow(outputVoltage, 2) / loadResistance;

            // Thermal noise power (kTBF model, simplified)
            double k = 1.38e-23; // Boltzmann constant (J/K)
            double T = 300;      // Temperature (K, ~27°C)
            double B = filterCutoff; // Bandwidth (Hz, using filter cutoff)
            double thermalNoisePower = k * T * B;

            // Additional noise from rectifier (simplified model)
            double rectifierNoiseFactor = rectifierType == "Full-Wave" ? 1.2 : 1.5;
            double totalNoisePower = thermalNoisePower * rectifierNoiseFactor;

            // SNR in dB
            double snr = 10 * Math.Log10(signalPower / totalNoisePower);
            return Math.Round(snr, 3);
        }

        public double CalculateTHD()
        {
            // Simplified THD model: harmonic distortion from rectifier and load
            double fundamentalAmplitude = outputVoltage;

            // Harmonic amplitudes (2nd and 3rd harmonics, simplified)
            double harmonicFactor = rectifierType == "Full-Wave" ? 0.05 : 0.1; // Full-wave: less distortion
            double secondHarmonic = fundamentalAmplitude * harmonicFactor;
            double thirdHarmonic = fundamentalAmplitude * harmonicFactor * 0.5;

            // Total harmonic power
            double harmonicPower = Math.Pow(secondHarmonic, 2) + Math.Pow(thirdHarmonic, 2);
            double fundamentalPower = Math.Pow(fundamentalAmplitude, 2);

            // THD as percentage
            double thd = Math.Sqrt(harmonicPower) / fundamentalPower * 100;
            return Math.Round(thd, 3);
        }

        public double CalculateEMI()
        {
            // Simplified conducted EMI model (dBµV) from PWM switching
            double baseEMI = 100; // Base EMI level (dBµV) at 1 kHz PWM
            double pwmEffect = 20 * Math.Log10(pwmCarrierFreq / 1000); // EMI scales with PWM freq
            double filterAttenuation = 10 * Math.Log10(1 + filterInductance * filterCapacitance * 1e-9 * Math.Pow(2 * Math.PI * filterCutoff, 2));

            // Total EMI (dBµV)
            double emi = baseEMI + pwmEffect - filterAttenuation;
            return Math.Round(emi, 3);
        }
    }
}