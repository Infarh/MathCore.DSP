using System;
using System.Linq;
using MathCore.DSP.Fourier;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using static System.Math;

namespace MathCore.DSP.Tests.Fourier
{
    [TestClass]
    public sealed class FFT_Test
    {
        /* ------------------------------------------------------------------------------------------ */

        //public TestContext TestContext { get; set; }

        #region Additional test attributes

        // 
        //You can use the following additional attributes as you write your tests:
        //

        ////Use ClassInitialize to run code before running the first test in the class
        //[ClassInitialize] public static void MyClassInitialize(TestContext context) { }

        //Use ClassCleanup to run code after all tests in a class have run
        //[ClassCleanup] public static void MyClassCleanup() { }

        //Use TestInitialize to run code before running each test
        //[TestInitialize] public void MyTestInitialize() { }

        //Use TestCleanup to run code after each test has run
        //[TestCleanup] public void MyTestCleanup() { }

        #endregion

        /* ------------------------------------------------------------------------------------------ */

        //private static void Clear(double[] array, double threshold = 1e-10)
        //{
        //    for (var i = 0; i < array.Length; i++)
        //        if (array[i] < threshold) array[i] = 0;
        //}

        private static void Clear(Complex[] array, double threshold = 1e-10)
        {
            for (var i = 0; i < array.Length; i++)
                if (array[i].Abs < threshold) array[i] = 0;
        }

        private static double TestFunction(double x) => 1
            + 2 * Cos(2 * PI * x * 1)
            + 8 * Cos(2 * PI * x * 4)
            + 12 * Cos(2 * PI * x * 6);

        [TestMethod]
        public void TransformTest()
        {
            Func<double, double> f = TestFunction;
            const double fd = 16;
            const double dt = 1 / fd;
            const int samples_count = 16;
            var signal = f.Sampling(0, dt, samples_count).ToArray();

            var spectrum = signal.FourierTransform();
            Clear(spectrum);
            var spectrum_abs = spectrum.ToAbs() ?? throw new AssertInconclusiveException();
            var spectrum_arg = spectrum.ToArgDeg() ?? throw new AssertInconclusiveException();
            Assert.That.Value(spectrum_abs.Length).IsEqual(samples_count);
            Assert.That.Value(spectrum_abs[0]).IsEqual(1, 6.45e-15);
            Assert.That.Value(spectrum_abs[1]).IsEqual(1, 5.56e-15);
            Assert.That.Value(spectrum_abs[4]).IsEqual(4, 1.79e-15);
            Assert.That.Value(spectrum_abs[6]).IsEqual(6, 3.56e-15);
            Assert.That.Value(spectrum_abs[samples_count - 6]).IsEqual(6, 4.45e-15);
            Assert.That.Value(spectrum_abs[samples_count-4]).IsEqual(4, 2.67e-15);
            Assert.That.Value(spectrum_abs[samples_count-1]).IsEqual(1, 5.23e-15);
            Assert.That.Value(spectrum_abs.Sum()).IsEqual(23, 3.56e-15);
            spectrum_arg.Foreach(v => Assert.That.Value(v).IsEqual(0, 1.47e-13));

            var sp = fft.FFT(signal);
            Clear(sp);
            sp.ToAbsArgDeg(out var sp_abs, out var sp_arg);
            Assert.That.Value(sp_abs.Length).IsEqual(samples_count);
            Assert.That.Value(sp_abs[0]).IsEqual(1, 6.45e-15);
            Assert.That.Value(sp_abs[1]).IsEqual(1, 5.56e-15);
            Assert.That.Value(sp_abs[4]).IsEqual(4, 1.79e-15);
            Assert.That.Value(sp_abs[6]).IsEqual(6, 3.56e-15);
            Assert.That.Value(sp_abs[samples_count - 6]).IsEqual(6, 4.45e-15);
            Assert.That.Value(sp_abs[samples_count - 4]).IsEqual(4, 2.67e-15);
            Assert.That.Value(sp_abs[samples_count - 1]).IsEqual(1, 5.23e-15);
            Assert.That.Value(sp_abs.Sum()).IsEqual(23, 3.56e-15);
            sp_arg.Foreach(v => Assert.That.Value(v).IsEqual(0, 1.44e-13));
        }

    }
}
