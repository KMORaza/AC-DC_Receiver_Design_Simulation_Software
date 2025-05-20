using System;

namespace ACDCSimulator
{
    public class ImpedanceMatcher
    {
        private double inputImpedance;
        private double outputImpedance;
        private double sourceImpedance;
        private double loadImpedance;

        public ImpedanceMatcher(double inputImpedance, double outputImpedance, double sourceImpedance, double loadImpedance)
        {
            this.inputImpedance = inputImpedance; // in ohms
            this.outputImpedance = outputImpedance; // in ohms
            this.sourceImpedance = sourceImpedance; // in ohms (e.g., 50 ohms)
            this.loadImpedance = loadImpedance; // in ohms (e.g., load resistance)
        }

        public double CalculateReflectionCoefficient()
        {
            // Reflection coefficient (Gamma) = |(Z_L - Z_S)/(Z_L + Z_S)|
            double numerator = Math.Abs(inputImpedance - sourceImpedance);
            double denominator = inputImpedance + sourceImpedance;
            double gamma = denominator != 0 ? numerator / denominator : 1.0;
            return Math.Round(gamma, 3);
        }

        public double CalculateMismatchLoss()
        {
            // Mismatch loss (dB) = -10 * log10(1 - Gamma^2)
            double gamma = CalculateReflectionCoefficient();
            double powerRatio = 1 - Math.Pow(gamma, 2);
            double mismatchLoss = powerRatio > 0 ? -10 * Math.Log10(powerRatio) : double.PositiveInfinity;
            return Math.Round(mismatchLoss, 3);
        }
    }
}