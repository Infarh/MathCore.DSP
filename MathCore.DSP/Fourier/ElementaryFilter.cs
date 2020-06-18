using System;
// ReSharper disable UnusedType.Global
// ReSharper disable NotAccessedField.Local
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Fourier
{
    /// <summary>Элементарный фильтр Фурье</summary>
    public class ElementaryFilter
    {
        /// <summary>Действительное значение</summary>
        private double _ReValue;
        /// <summary>Мнимое значение</summary>
        private double _ImValue;
        /// <summary>Размер выборки</summary>
        private readonly int _N;
        /// <summary>Количество спектральных компонент в спектре</summary>
#pragma warning disable IDE0052 // Удалить непрочитанные закрытые члены
        private readonly int _M;
#pragma warning restore IDE0052 // Удалить непрочитанные закрытые члены
        /// <summary>Дискрет фазы</summary>
        private readonly double _dPhi;
        /// <summary>Индекс отсчёта</summary>
        private int _SamplesIndex;

        /// <summary>Значение фильтра</summary>
        public Complex Value => new Complex(_ReValue / _N, _ImValue / _N);

        /// <summary>Инициализация нового элементарного фильтра преобразования Фурье</summary>
        /// <param name="N">Размер выборки</param>
        /// <param name="M">Размер спектра</param>
        public ElementaryFilter(int N, int M)
        {
            _M = M;
            _N = N;
            _dPhi = M * Consts.pi2 / N;
        }

        /// <summary>Инициализация фильтра (сброс состояния)</summary>
        public void Initialize()
        {
            _ReValue = 0;
            _ImValue = 0;
            _SamplesIndex = 0;
        }

        /// <summary>Обработка очередного значения</summary>
        /// <param name="value">Значение обрабатываемого отсчёта</param>
        public Complex Process(double value)
        {
            var arg = _dPhi * _SamplesIndex++;
            _ReValue += Math.Cos(arg) * value;
            _ImValue += Math.Sin(arg) * value;
            return Value;
        }
    }
}
