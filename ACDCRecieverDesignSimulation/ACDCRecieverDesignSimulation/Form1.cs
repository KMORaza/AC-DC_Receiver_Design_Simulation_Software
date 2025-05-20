using System;
using System.Drawing;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using System.Linq;

namespace ACDCSimulator
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
            // Button hover effects
            btnSimulate.MouseEnter += (s, e) => btnSimulate.BackColor = Color.FromArgb(192, 192, 192);
            btnSimulate.MouseLeave += (s, e) => btnSimulate.BackColor = Color.FromArgb(160, 160, 160);
            // Populate dropdowns
            comboRectifier.Items.AddRange(new object[] { RectifierType.HalfWave, RectifierType.FullWave, RectifierType.Bridge });
            comboFilter.Items.AddRange(new object[] { FilterType.None, FilterType.Capacitive, FilterType.Inductive, FilterType.Active });
            comboRegulator.Items.AddRange(new object[] { RegulatorType.None, RegulatorType.Linear, RegulatorType.Switching });
            comboSwitchingDevice.Items.AddRange(new object[] { SwitchingDeviceType.MOSFET, SwitchingDeviceType.IGBT, SwitchingDeviceType.GaN, SwitchingDeviceType.SiC });
            comboRectifier.SelectedIndex = 2;
            comboFilter.SelectedIndex = 0;
            comboRegulator.SelectedIndex = 0;
            comboSwitchingDevice.SelectedIndex = 0;
            // Clear charts to prevent initialization issues
            chartBode.Series.Clear();
            chartNyquist.Series.Clear();
            chartRootLocus.Series.Clear();
            chartWaveform.Series.Clear();
            // Configure chartBode to prevent logarithmic scale errors
            chartBode.ChartAreas[0].AxisX.IsLogarithmic = false; // Disable temporarily to avoid errors
            chartBode.ChartAreas[0].AxisX.Minimum = 0.1;
            chartBode.ChartAreas[0].AxisX.Maximum = 100000; // Reasonable max (100 kHz)
            chartBode.Visible = false; // Hide chart until valid data is available
            // Handle tab selection to show chart only with valid data
            tabControl.SelectedIndexChanged += (s, e) =>
            {
                if (tabControl.SelectedTab == tabAnalysis)
                {
                    chartBode.Visible = chartBode.Series.Count > 0;
                }
            };
        }

        private async void btnSimulate_Click(object sender, EventArgs e)
        {
            try
            {
                // Validate inputs
                if (!double.TryParse(txtInputVoltage.Text, out double inputVoltage) || inputVoltage <= 0)
                {
                    MessageBox.Show("Please enter a valid positive input AC voltage.", "Input Error");
                    return;
                }
                if (!double.TryParse(txtLoadResistance.Text, out double loadResistance) || loadResistance <= 0)
                {
                    MessageBox.Show("Please enter a valid positive load resistance.", "Input Error");
                    return;
                }
                if (!double.TryParse(txtTurnsRatio.Text, out double turnsRatio) || turnsRatio <= 0)
                {
                    MessageBox.Show("Please enter a valid positive turns ratio.", "Input Error");
                    return;
                }
                double? filterCapacitance = null, filterInductance = null, filterCutoff = null, regulatorVoltage = null;
                if (comboFilter.SelectedItem?.ToString() == FilterType.Capacitive.ToString())
                {
                    if (!double.TryParse(txtFilterCapacitance.Text, out double cap) || cap <= 0)
                    {
                        MessageBox.Show("Please enter a valid positive capacitance.", "Input Error");
                        return;
                    }
                    filterCapacitance = cap * 1e-6;
                }
                else if (comboFilter.SelectedItem?.ToString() == FilterType.Inductive.ToString())
                {
                    if (!double.TryParse(txtFilterInductance.Text, out double ind) || ind <= 0)
                    {
                        MessageBox.Show("Please enter a valid positive inductance.", "Input Error");
                        return;
                    }
                    filterInductance = ind * 1e-3;
                }
                else if (comboFilter.SelectedItem?.ToString() == FilterType.Active.ToString())
                {
                    if (!double.TryParse(txtFilterCutoff.Text, out double freq) || freq <= 0)
                    {
                        MessageBox.Show("Please enter a valid positive cutoff frequency.", "Input Error");
                        return;
                    }
                    filterCutoff = freq;
                }
                if (comboRegulator.SelectedItem?.ToString() == RegulatorType.Linear.ToString() ||
                    comboRegulator.SelectedItem?.ToString() == RegulatorType.Switching.ToString())
                {
                    if (!double.TryParse(txtRegulatorVoltage.Text, out double volt) || volt <= 0)
                    {
                        MessageBox.Show("Please enter a valid positive regulator voltage.", "Input Error");
                        return;
                    }
                    regulatorVoltage = volt;
                }
                if (!double.TryParse(txtPIDKp.Text, out double pidKp) || pidKp < 0)
                {
                    MessageBox.Show("Please enter a valid non-negative PID Kp.", "Input Error");
                    return;
                }
                if (!double.TryParse(txtPIDKi.Text, out double pidKi) || pidKi < 0)
                {
                    MessageBox.Show("Please enter a valid non-negative PID Ki.", "Input Error");
                    return;
                }
                if (!double.TryParse(txtPIDKd.Text, out double pidKd) || pidKd < 0)
                {
                    MessageBox.Show("Please enter a valid non-negative PID Kd.", "Input Error");
                    return;
                }
                if (!double.TryParse(txtPLLRefFreq.Text, out double pllRefFreq) || pllRefFreq <= 0)
                {
                    MessageBox.Show("Please enter a valid positive PLL reference frequency.", "Input Error");
                    return;
                }
                if (!double.TryParse(txtPWMCarrierFreq.Text, out double pwmCarrierFreq) || pwmCarrierFreq <= 0)
                {
                    MessageBox.Show("Please enter a valid positive PWM carrier frequency.", "Input Error");
                    return;
                }
                double[] digitalCoefficients;
                try
                {
                    digitalCoefficients = txtDigitalCoeff.Text.Split(new[] { ',' }, StringSplitOptions.RemoveEmptyEntries)
                                                         .Select(s => double.Parse(s.Trim()))
                                                         .ToArray();
                    if (digitalCoefficients.Length == 0)
                    {
                        MessageBox.Show("Please enter at least one valid digital coefficient.", "Input Error");
                        return;
                    }
                }
                catch (FormatException)
                {
                    MessageBox.Show("Digital coefficients must be valid numbers.", "Input Error");
                    return;
                }

                // Disable button and show progress
                btnSimulate.Enabled = false;
                lblStatus.Text = "Simulating...";
                lblStatus.ForeColor = Color.White;

                // Create circuit parameters
                var parameters = new CircuitParameters
                {
                    Rectifier = (RectifierType)comboRectifier.SelectedItem,
                    Filter = (FilterType)comboFilter.SelectedItem,
                    Regulator = (RegulatorType)comboRegulator.SelectedItem,
                    SwitchingDevice = (SwitchingDeviceType)comboSwitchingDevice.SelectedItem,
                    TransformerTurnsRatio = turnsRatio,
                    FilterCapacitance = filterCapacitance ?? 100e-6,
                    FilterInductance = filterInductance ?? 10e-3,
                    FilterCutoffFrequency = filterCutoff ?? 100,
                    RegulatorVoltage = regulatorVoltage ?? 5.0,
                    PID_Kp = pidKp,
                    PID_Ki = pidKi,
                    PID_Kd = pidKd,
                    PLL_ReferenceFrequency = pllRefFreq,
                    PWM_CarrierFrequency = pwmCarrierFreq,
                    DigitalCoefficients = digitalCoefficients
                };

                // Run simulation
                SimulationResult result = await Task.Run(() =>
                {
                    var simulator = new CircuitSimulator(inputVoltage, loadResistance, parameters);
                    return simulator.RunSimulation();
                });

                // Calculate additional metrics
                // 1. Compute Input and Output Impedance based on configuration
                double inputImpedance = loadResistance; // Base impedance
                double outputImpedance = loadResistance; // Base impedance
                double sourceImpedance = 50; // Typical source impedance

                // Adjust input impedance based on rectifier and transformer
                double rectifierFactor = 1.0;
                if (parameters.Rectifier == RectifierType.HalfWave)
                    rectifierFactor = 0.5;
                else if (parameters.Rectifier == RectifierType.FullWave)
                    rectifierFactor = 0.8;
                else if (parameters.Rectifier == RectifierType.Bridge)
                    rectifierFactor = 1.0;

                // Reflect load impedance to primary side of transformer
                inputImpedance = loadResistance / (parameters.TransformerTurnsRatio * parameters.TransformerTurnsRatio) * rectifierFactor;

                // Adjust output impedance based on filter and regulator
                double filterFactor = 1.0;
                if (parameters.Filter == FilterType.Capacitive)
                {
                    // Simplified: Capacitor reduces impedance at higher frequencies
                    double freq = 60; // Default frequency
                    double xc = 1 / (2 * Math.PI * freq * parameters.FilterCapacitance);
                    filterFactor = 1 / (1 + loadResistance / xc);
                }
                else if (parameters.Filter == FilterType.Inductive)
                {
                    // Simplified: Inductor increases impedance at higher frequencies
                    double freq = 60;
                    double xl = 2 * Math.PI * freq * parameters.FilterInductance;
                    filterFactor = 1 + xl / loadResistance;
                }
                else if (parameters.Filter == FilterType.Active)
                {
                    // Active filter: Assume it maintains impedance close to load
                    filterFactor = 1.0;
                }

                double regulatorFactor = 1.0;
                if (parameters.Regulator == RegulatorType.Linear)
                    regulatorFactor = 0.9; // Linear regulators slightly reduce impedance
                else if (parameters.Regulator == RegulatorType.Switching)
                    regulatorFactor = 0.7; // Switching regulators reduce impedance more

                outputImpedance = loadResistance * filterFactor * regulatorFactor;

                // 2. Impedance Matching (Reflection Coefficient, Mismatch Loss)
                var impedanceMatcher = new ImpedanceMatcher(
                    inputImpedance: inputImpedance,
                    outputImpedance: outputImpedance,
                    sourceImpedance: sourceImpedance,
                    loadImpedance: loadResistance
                );
                double reflectionCoefficient = impedanceMatcher.CalculateReflectionCoefficient();
                double mismatchLoss = impedanceMatcher.CalculateMismatchLoss();

                // 3. Noise and Distortion (SNR, EMI, THD)
                var noiseAnalyzer = new NoiseAndDistortionAnalyzer(
                    inputVoltage: inputVoltage,
                    loadResistance: loadResistance,
                    rectifierType: comboRectifier.SelectedItem.ToString(),
                    filterCapacitance: filterCapacitance ?? 100e-6,
                    filterInductance: filterInductance ?? 10e-3,
                    filterCutoff: filterCutoff ?? 100,
                    pwmCarrierFreq: pwmCarrierFreq,
                    outputVoltage: result.SteadyStateDcVoltage
                );
                double snr = noiseAnalyzer.CalculateSNR();
                double emi = noiseAnalyzer.CalculateEMI();
                double thd = noiseAnalyzer.CalculateTHD();

                // 4. Power Factor
                var powerFactorCorrector = new PowerFactorCorrector(
                    inputVoltage: inputVoltage,
                    loadResistance: loadResistance,
                    rectifierType: comboRectifier.SelectedItem.ToString(),
                    frequency: 60 // Default frequency as per CircuitSimulator.cs
                );
                double powerFactor = powerFactorCorrector.ApplyPowerFactorCorrection();

                // Update results
                lblOutputVoltage.Text = $"{result.SteadyStateDcVoltage:F2} V";
                lblRippleVoltage.Text = $"{result.RippleVoltage:F2} V";
                lblPowerDissipation.Text = $"{result.PowerDissipation:F2} W";
                lblDominantHarmonic.Text = $"{result.DominantHarmonicFrequency:F2} Hz, {result.DominantHarmonicAmplitude:F2} V";
                lblTotalHarmonicDistortion.Text = $"{thd:F2} %";
                lblPowerFactor.Text = $"{powerFactor:F3}";
                lblSignalToNoiseRatio.Text = $"{snr:F2} dB";
                lblEMI.Text = $"{emi:F2} dBµV";
                lblInputImpedance.Text = $"{inputImpedance:F2} Ω";
                lblOutputImpedance.Text = $"{outputImpedance:F2} Ω";
                lblReflectionCoefficient.Text = $"{reflectionCoefficient:F3}";
                lblMismatchLoss.Text = $"{mismatchLoss:F2} dB";

                // Plot transient waveform
                chartWaveform.Series.Clear();
                var series = new Series("Output Voltage")
                {
                    ChartType = SeriesChartType.Line,
                    Color = Color.Black
                };
                for (int i = 0; i < result.TransientTimePoints.Length; i++)
                {
                    series.Points.AddXY(result.TransientTimePoints[i] * 1000, result.TransientVoltagePoints[i]);
                }
                chartWaveform.Series.Add(series);

                // Plot Bode plot
                chartBode.Series.Clear();
                var magSeries = new Series("Magnitude (dB)")
                {
                    ChartType = SeriesChartType.Line,
                    Color = Color.Black
                };
                var phaseSeries = new Series("Phase (deg)")
                {
                    ChartType = SeriesChartType.Line,
                    Color = Color.DarkGray
                };
                for (int i = 0; i < result.BodeFrequencyPoints.Length; i++)
                {
                    if (result.BodeFrequencyPoints[i] > 0) // Only plot positive frequencies
                    {
                        magSeries.Points.AddXY(result.BodeFrequencyPoints[i], result.BodeMagnitudePoints[i]);
                        phaseSeries.Points.AddXY(result.BodeFrequencyPoints[i], result.BodePhasePoints[i]);
                    }
                }
                if (magSeries.Points.Count > 0) // Only add series if valid data exists
                {
                    chartBode.Series.Add(magSeries);
                    chartBode.Series.Add(phaseSeries);
                    chartBode.ChartAreas[0].AxisX.IsLogarithmic = true; // Re-enable logarithmic scale
                    chartBode.Visible = true; // Show chart after valid data
                }
                else
                {
                    MessageBox.Show("No valid frequency points for Bode plot.", "Plotting Warning");
                    chartBode.Visible = false;
                }

                // Plot Nyquist plot
                chartNyquist.Series.Clear();
                var nyquistSeries = new Series("Nyquist")
                {
                    ChartType = SeriesChartType.Point,
                    Color = Color.Black
                };
                for (int i = 0; i < result.NyquistRealPoints.Length; i++)
                {
                    nyquistSeries.Points.AddXY(result.NyquistRealPoints[i], result.NyquistImagPoints[i]);
                }
                chartNyquist.Series.Add(nyquistSeries);

                // Plot Root Locus
                chartRootLocus.Series.Clear();
                var rootSeries = new Series("Root Locus")
                {
                    ChartType = SeriesChartType.Point,
                    Color = Color.Black
                };
                foreach (var root in result.RootLocusPoints)
                {
                    rootSeries.Points.AddXY(root.Real, root.Imaginary);
                }
                chartRootLocus.Series.Add(rootSeries);

                lblStatus.Text = "Simulation Complete";
                lblStatus.ForeColor = Color.LimeGreen;
            }
            catch (Exception ex)
            {
                MessageBox.Show($"An error occurred: {ex.Message}", "Simulation Error");
                lblStatus.Text = "Simulation Failed";
                lblStatus.ForeColor = Color.Red;
            }
            finally
            {
                btnSimulate.Enabled = true;
            }
        }
    }
}