using System;
using System.Collections.Generic;
using System.Text.RegularExpressions;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Shapes;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Legends;
using OxyPlot.Series;
using OxyPlot.Wpf;

namespace AsteroidOrbitSimulation
{
    public partial class MainWindow : Window
    {
        private PlotModel plotModel;
        private List<DataPoint> originalOrbitPoints = new List<DataPoint>();
        private List<DataPoint> currentTrajectoryPoints = new List<DataPoint>();
        private Simulation simulation;
        private bool isSimulating;
        private Ellipse asteroidEllipse;

        public MainWindow()
        {
            InitializeComponent();
            InitializePlot();
            simulation = new Simulation();
        }

        private void InitializePlot()
        {
            plotModel = new PlotModel
            {
                Title = "Сравнение орбит",
                PlotAreaBorderColor = OxyColors.Black
            };

            plotModel.Legends.Add(new Legend
            {
                LegendPosition = LegendPosition.RightTop,
                LegendOrientation = LegendOrientation.Vertical,
                FontSize = 12
            });

            plotModel.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Bottom,
                Title = "X (а.е.)",
                MajorGridlineStyle = LineStyle.Solid,
                MinorGridlineStyle = LineStyle.Dot,
                Minimum = -2,
                Maximum = 2
            });

            plotModel.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Left,
                Title = "Y (а.е.)",
                MajorGridlineStyle = LineStyle.Solid,
                MinorGridlineStyle = LineStyle.Dot,
                Minimum = -2,
                Maximum = 2
            });

            var gridLines = new LineSeries
            {
                Color = OxyColors.LightGray,
                StrokeThickness = 0.5,
                LineStyle = LineStyle.Dot
            };

            for (double x = -5; x <= 5; x++)
            {
                gridLines.Points.Add(new DataPoint(x, -5));
                gridLines.Points.Add(new DataPoint(x, 5));
            }

            for (double y = -5; y <= 5; y++)
            {
                gridLines.Points.Add(new DataPoint(-5, y));
                gridLines.Points.Add(new DataPoint(5, y));
            }

            plotModel.Series.Add(gridLines);

            var originalOrbit = new LineSeries
            {
                Title = "Исходная орбита",
                Color = OxyColors.Navy,
                LineStyle = LineStyle.Dash,
                StrokeThickness = 2
            };

            var currentTrajectory = new LineSeries
            {
                Title = "Текущая траектория",
                Color = OxyColors.DarkRed,
                StrokeThickness = 2
            };

            plotModel.Series.Add(originalOrbit);
            plotModel.Series.Add(currentTrajectory);

            PlotView.Model = plotModel;
        }

        private void StartSimulation_Click(object sender, RoutedEventArgs e)
        {
            if (isSimulating) return;

            if (!TryParseInput(EjectionTimeInput.Text, out double ejectionTime) ||
                !TryParseInput(EjectionForceInput.Text, out double ejectionForce) ||
                !TryParseInput(EjectionAngleInput.Text, out double ejectionAngle))
            {
                MessageBox.Show("Проверьте корректность введенных параметров", "Ошибка ввода", MessageBoxButton.OK, MessageBoxImage.Error);
                return;
            }

            isSimulating = true;

            originalOrbitPoints.Clear();
            currentTrajectoryPoints.Clear();
            CanvasTrajectory.Children.Clear();

            simulation.Initialize(
                initialX: 1.0,
                initialY: 0.0,
                initialVx: 0.0,
                initialVy: 2 * Math.PI,
                ejectionTime: ejectionTime,
                ejectionForce: ejectionForce,
                ejectionAngle: ejectionAngle
            );

            originalOrbitPoints = GenerateOrbitPoints(
                simulation.InitialSemiMajorAxis,
                simulation.InitialEccentricity
            );

            (plotModel.Series[0] as LineSeries).Points.Clear();
            (plotModel.Series[0] as LineSeries).Points.AddRange(originalOrbitPoints);
            (plotModel.Series[1] as LineSeries).Points.Clear();

            var sun = new Ellipse
            {
                Width = 20,
                Height = 20,
                Fill = Brushes.Yellow
            };

            Canvas.SetLeft(sun, 500 - 10);
            Canvas.SetTop(sun, 300 - 10);
            CanvasTrajectory.Children.Add(sun);

            asteroidEllipse = new Ellipse
            {
                Width = 5,
                Height = 5,
                Fill = Brushes.Red
            };

            CanvasTrajectory.Children.Add(asteroidEllipse);

            CompositionTarget.Rendering += UpdateSimulation;
        }

        private bool TryParseInput(string input, out double result)
        {
            return double.TryParse(input.Replace(',', '.'),
                                  System.Globalization.NumberStyles.Any,
                                  System.Globalization.CultureInfo.InvariantCulture,
                                  out result);
        }

        private void NumberValidationTextBox(object sender, TextCompositionEventArgs e)
        {
            var regex = new Regex(@"^[0-9.,]$");
            e.Handled = !regex.IsMatch(e.Text);
        }

        private List<DataPoint> GenerateOrbitPoints(double semiMajorAxis, double eccentricity)
        {
            var points = new List<DataPoint>();
            double theta = 0;

            while (theta < 2 * Math.PI)
            {
                double r = (semiMajorAxis * (1 - Math.Pow(eccentricity, 2))) /
                           (1 + eccentricity * Math.Cos(theta));
                double x = r * Math.Cos(theta);
                double y = r * Math.Sin(theta);
                points.Add(new DataPoint(x, y));
                theta += Math.PI / 64;
            }

            return points;
        }

        private void UpdateSimulation(object sender, EventArgs e)
        {
            if (!isSimulating) return;

            simulation.Update(0.01);
            currentTrajectoryPoints.Add(new DataPoint(simulation.X, simulation.Y));

            var currentSeries = plotModel.Series[1] as LineSeries;
            currentSeries.Points.Clear();
            currentSeries.Points.AddRange(currentTrajectoryPoints);

            Canvas.SetLeft(asteroidEllipse, 500 + simulation.X * 100 - 2.5);
            Canvas.SetTop(asteroidEllipse, 300 - simulation.Y * 100 - 2.5);

            OrbitParams.Text = $"Начальная a: {simulation.InitialSemiMajorAxis:F2} а.е.\n" +
                               $"Текущая a: {simulation.CurrentSemiMajorAxis:F2} а.е.\n" +
                               $"Начальный e: {simulation.InitialEccentricity:F2}\n" +
                               $"Текущий e: {simulation.CurrentEccentricity:F2}\n" +
                               $"Скорость: {Math.Sqrt(simulation.Vx * simulation.Vx + simulation.Vy * simulation.Vy):F2} а.е./год\n" +
                               $"Расстояние: {Math.Sqrt(simulation.X * simulation.X + simulation.Y * simulation.Y):F2} а.е.\n" +
                               $"Время: {simulation.Time:F1} лет";

            plotModel.InvalidatePlot(true);

            if (simulation.Time > 20) StopSimulation();
        }

        private void StopSimulation()
        {
            isSimulating = false;
            CompositionTarget.Rendering -= UpdateSimulation;
        }

        private void StopButton_Click(object sender, RoutedEventArgs e) => StopSimulation();
    }

    public class Simulation
    {
        private const double G = 4 * Math.PI * Math.PI;
        private const double SunMass = 1.0;

        public double X { get; private set; }
        public double Y { get; private set; }
        public double Vx { get; private set; }
        public double Vy { get; private set; }
        public double Time { get; private set; }
        public double InitialSemiMajorAxis { get; private set; }
        public double InitialEccentricity { get; private set; }

        private double ejectionTime;
        private double ejectionForce;
        private double ejectionAngle;

        public void Initialize(double initialX, double initialY, double initialVx, double initialVy,
            double ejectionTime, double ejectionForce, double ejectionAngle)
        {
            X = initialX;
            Y = initialY;
            Vx = initialVx;
            Vy = initialVy;
            Time = 0;
            this.ejectionTime = ejectionTime;
            this.ejectionForce = ejectionForce;
            this.ejectionAngle = ejectionAngle * Math.PI / 180;

            CalculateInitialOrbitParameters();
        }

        private void CalculateInitialOrbitParameters()
        {
            double r = Math.Sqrt(X * X + Y * Y);
            double vSquared = Vx * Vx + Vy * Vy;
            double energy = 0.5 * vSquared - G * SunMass / r;
            InitialSemiMajorAxis = -G * SunMass / (2 * energy);
            double h = X * Vy - Y * Vx;
            InitialEccentricity = Math.Sqrt(1 + (2 * energy * h * h) / (G * G * SunMass * SunMass));
        }

        public void Update(double dt)
        {
            double r = Math.Sqrt(X * X + Y * Y);
            double ax = -G * SunMass * X / Math.Pow(r, 3);
            double ay = -G * SunMass * Y / Math.Pow(r, 3);

            if (Math.Abs(Time - ejectionTime) < dt / 2)
            {
                ax += ejectionForce * Math.Cos(ejectionAngle);
                ay += ejectionForce * Math.Sin(ejectionAngle);
            }

            Vx += ax * dt;
            Vy += ay * dt;
            X += Vx * dt;
            Y += Vy * dt;
            Time += dt;
        }

        public double CurrentSemiMajorAxis
        {
            get
            {
                double r = Math.Sqrt(X * X + Y * Y);
                double vSquared = Vx * Vx + Vy * Vy;
                double energy = 0.5 * vSquared - G * SunMass / r;
                return energy < 0 ? -G * SunMass / (2 * energy) : double.PositiveInfinity;
            }
        }

        public double CurrentEccentricity
        {
            get
            {
                double r = Math.Sqrt(X * X + Y * Y);
                double vSquared = Vx * Vx + Vy * Vy;
                double energy = 0.5 * vSquared - G * SunMass / r;

                if (energy >= 0) return 1.0;

                double h = X * Vy - Y * Vx;
                return Math.Sqrt(1 + (2 * energy * h * h) / (G * G * SunMass * SunMass));
            }
        }
    }
}