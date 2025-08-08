using MathCore.DSP.Samples;

namespace MathCore.DSP.Tests.Samples;

[TestClass]
public class SampleSI16CalculatorTests
{
    /// <summary>Вспомогательный класс для доступа к защищённому методу GetIndex</summary>
    private class SampleSI16CalculatorWrapper : SampleSI16Calculator
    {
        /// <summary>Публичный метод для тестирования GetIndex</summary>
        public int GetIndexPublic(int I, int Q) => GetIndex(I, Q); // Получить индекс
    }

    /// <summary>Проверяет уникальность и корректность индексов для всех пар I, Q в диапазоне [0, 127]</summary>
    [TestMethod]
    public void GetIndex_AllCombinations_UniqueAndSymmetric()
    {
        var wrapper = new SampleSI16CalculatorWrapper();
        var indices = new HashSet<int>();
        for (var i = 0; i <= 128; i++)
            for (var q = 0; q <= 128; q++)
            {
                var index1 = wrapper.GetIndexPublic(i, q);
                var index2 = wrapper.GetIndexPublic(q, i);

                Assert.AreEqual(index1, index2); // Проверка симметрии
                Assert.IsTrue(index1 is >= 0 and <= 8384); // Проверка диапазона

                indices.Add(index1);
            }

        Assert.AreEqual(8385, indices.Count); // Проверка уникальности
    }
}
