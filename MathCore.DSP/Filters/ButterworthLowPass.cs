using static System.Math;

using static MathCore.Polynom.Array;
using static MathCore.SpecialFunctions;

namespace MathCore.DSP.Filters;

/// <summary>Низкочастотный фильтр Баттерворта</summary>
/// <remarks>
/// Фильтр Баттерворта - это фильтр с максимально плоской амплитудной характеристикой в полосе пропускания.
/// Реализует цифровой фильтр на основе аналогового прототипа с использованием билинейного преобразования.
/// Фильтр характеризуется хорошей прямоугольностью АЧХ и быстрой переходной характеристикой.
/// </remarks>
public class ButterworthLowPass : ButterworthFilter
{
    /// <summary>Вычисляет минимальный порядок фильтра по требуемым параметрам частотной характеристики</summary>
    /// <param name="dt">Период дискретизации (обратное значение частоты дискретизации), секунды</param>
    /// <param name="fp">Граничная частота полосы пропускания, Гц</param>
    /// <param name="fs">Граничная частота полосы заграждения (подавления), Гц</param>
    /// <param name="Gp">Коэффициент передачи в полосе пропускания, по умолчанию 0.891250938 (эквивалент -1 дБ)</param>
    /// <param name="Gs">Коэффициент передачи в полосе заграждения, по умолчанию 0.01 (эквивалент -40 дБ)</param>
    /// <returns>Минимальный порядок фильтра для выполнения спецификации</returns>
    /// <remarks>
    /// Метод использует формулу расчета порядка фильтра Баттерворта на основе логарифмического отношения
    /// коэффициентов затухания и отношения частот. Результат округляется в большую сторону.
    /// </remarks>
    /// <example>
    /// <![CDATA[
    /// // Расчет порядка фильтра с частотой дискретизации 1000 Гц,
    /// // полосой пропускания до 100 Гц и полосой заграждения от 200 Гц
    /// double dt = 0.001; // 1000 Гц
    /// int order = ButterworthLowPass.GetOrder(dt, fp: 100, fs: 200);
    /// // Результат: order >= требуемому для выполнения спецификации
    /// ]]>
    /// </example>
    public static int GetOrder(double dt, double fp, double fs, double Gp = 0.891250938, double Gs = 0.01)
    {
        var kEps2 = (1 / (Gs * Gs) - 1) / (1 / (Gp * Gp) - 1);
        var pi_dt = PI * dt;
        var kW = Tan(fs * pi_dt) / Tan(fp * pi_dt);

        var N = (int)Ceiling(0.5 * Log(kEps2) / Log(kW));
        return N;
    }

    /// <summary>Вычисляет граничную частоту полосы заграждения по заданному порядку фильтра</summary>
    /// <param name="dt">Период дискретизации, секунды</param>
    /// <param name="fp">Граничная частота полосы пропускания, Гц</param>
    /// <param name="Order">Порядок фильтра</param>
    /// <param name="Gp">Коэффициент передачи в полосе пропускания, по умолчанию 0.891250938 (-1 дБ)</param>
    /// <param name="Gs">Коэффициент передачи в полосе заграждения, по умолчанию 0.01 (-40 дБ)</param>
    /// <returns>Граничная частота полосы заграждения в Гц</returns>
    /// <remarks>
    /// Используется для определения граничной частоты заграждения, когда известны порядок фильтра и требуемые
    /// характеристики затухания. Полезно при итеративном проектировании фильтра.
    /// </remarks>
    /// <example>
    /// <![CDATA[
    /// // Определение частоты заграждения для фильтра 4-го порядка
    /// double dt = 0.001;
    /// double fs = ButterworthLowPass.GetFrequencyStop(dt, fp: 100, Order: 4);
    /// ]]>
    /// </example>
    public static double GetFrequencyStop(double dt, double fp, int Order, double Gp = 0.891250938, double Gs = 0.01)
    {
        var kEps2 = (1 / (Gs * Gs) - 1) / (1 / (Gp * Gp) - 1);
        var pi_dt = PI * dt;

        var fs = Atan(Pow(kEps2, 0.5 / Order) * Tan(pi_dt * fp)) / pi_dt;
        return fs;
    }

    /// <summary>Вычисляет граничную частоту полосы пропускания по заданному порядку фильтра</summary>
    /// <param name="dt">Период дискретизации, секунды</param>
    /// <param name="Order">Порядок фильтра</param>
    /// <param name="fs">Граничная частота полосы заграждения, Гц</param>
    /// <param name="Gp">Коэффициент передачи в полосе пропускания, по умолчанию 0.891250938 (-1 дБ)</param>
    /// <param name="Gs">Коэффициент передачи в полосе заграждения, по умолчанию 0.01 (-40 дБ)</param>
    /// <returns>Граничная частота полосы пропускания в Гц</returns>
    /// <remarks>
    /// Используется для определения граничной частоты пропускания, когда известны порядок фильтра и требуемые
    /// характеристики затухания. Применяется при обратном расчете спецификации фильтра.
    /// </remarks>
    /// <example>
    /// <![CDATA[
    /// // Определение частоты пропускания для фильтра 4-го порядка
    /// double dt = 0.001;
    /// double fp = ButterworthLowPass.GetFrequencyPass(dt, Order: 4, fs: 200);
    /// ]]>
    /// </example>
    public static double GetFrequencyPass(double dt, int Order, double fs, double Gp = 0.891250938, double Gs = 0.01)
    {
        var kEps2 = (1 / (Gs * Gs) - 1) / (1 / (Gp * Gp) - 1);
        var pi_dt = PI * dt;

        var fp = Atan(Pow(kEps2, -0.5 / Order) * Tan(pi_dt * fs)) / pi_dt;
        return fp;
    }

    /// <summary>Получает набор полюсов нормированного фильтра Баттерворта</summary>
    /// <param name="Order">Порядок фильтра (количество полюсов)</param>
    /// <param name="Gp">Коэффициент передачи, по умолчанию 0.891250938 (-1 дБ)</param>
    /// <returns>Перечисление комплексных полюсов нормированного фильтра</returns>
    /// <remarks>
    /// Полюса вычисляются с использованием формулы полюсов Баттерворта в s-плоскости.
    /// Нормирование означает, что фильтр рассчитан на частоту 1 рад/с и требует дальнейшего
    /// преобразования для получения требуемой частоты среза.
    /// </remarks>
    /// <example>
    /// <![CDATA[
    /// // Получение полюсов фильтра 3-го порядка
    /// var poles = ButterworthLowPass.GetPoles(Order: 3);
    /// foreach(var pole in poles)
    /// {
    ///     Console.WriteLine($"Полюс: {pole}");
    /// }
    /// ]]>
    /// </example>
    public static IEnumerable<Complex> GetPoles(int Order, double Gp = 0.891250938) => GetNormPolesGp(Order, Gp);

    // https://ru.dsplib.org/content/filter_butter_ap/filter_butter_ap.html

    /// <summary>Вычисляет коэффициенты передаточной функции фильтра Баттерворта</summary>
    /// <param name="Spec">Спецификация фильтра с требуемыми параметрами частотной характеристики</param>
    /// <returns>Кортеж (A, B) где A - коэффициенты полинома знаменателя, B - коэффициенты полинома числителя передаточной функции</returns>
    /// <remarks>
    /// Метод выполняет следующие этапы синтеза фильтра:
    /// 1. Вычисляет порядок фильтра на основе спецификации
    /// 2. Получает нормированные полюса фильтра Баттерворта
    /// 3. Трансформирует полюса на требуемую частоту пропускания
    /// 4. Преобразует полюса из s-плоскости (аналоговый фильтр) в z-плоскость (цифровой фильтр) 
    ///    с использованием билинейного преобразования
    /// 5. Вычисляет нормирующий коэффициент
    /// 6. Формирует коэффициенты полиномов числителя и знаменателя
    /// 
    /// Коэффициенты A и B используются для рекуррентного разностного уравнения фильтра.
    /// </remarks>
    /// <example>
    /// <![CDATA[
    /// // Расчет коэффициентов фильтра Баттерворта
    /// var spec = new Specification(dt: 0.001, fp: 100, fs: 200, Gp: 0.891250938, Gs: 0.01);
    /// var (A, B) = ButterworthLowPass.GetPolynoms(spec);
    /// 
    /// // A содержит коэффициенты знаменателя (например, [1.0, -2.5, 2.1, -0.6])
    /// // B содержит коэффициенты числителя (например, [0.05, 0.15, 0.15, 0.05])
    /// ]]>
    /// </example>
    public static (double[] A, double[] B) GetPolynoms(Specification Spec)
    {
        // Порядок фильтра
        var N = (int)Ceiling(Log(Spec.kEps) / Log(Spec.kW));
        var poles = GetNormPoles(N, Spec.EpsP);

        // Масштабируем полюса на требуемую частоту пропускания
        var Wp = Spec.Wp;
        var translated_poles = TransformToLowPassW(poles, Wp);

        // Переходим из p-плоскости в z-плоскость
        var dt = Spec.dt;
        var z_poles = ToZArray(translated_poles, dt);

        // Вычисляем нормирующий множитель
        var g_norm = z_poles.Multiply(z => (1 - z) / 2).Re;

        var B = new double[N + 1];
        for (var i = 0; i < B.Length; i++)
            B[i] = g_norm * BinomialCoefficient(N, i);

        var A = GetCoefficientsInverted(z_poles).ToRe();

        return (A, B);
    }

    /// <summary>Инициализирует новый низкочастотный фильтр Баттерворта с указанными параметрами</summary>
    /// <param name="dt">Период дискретизации (обратное значение частоты дискретизации), секунды</param>
    /// <param name="fp">Граничная частота полосы пропускания, Гц. Должна быть меньше половины частоты дискретизации</param>
    /// <param name="fs">Граничная частота полосы заграждения, Гц. Должна быть больше частоты пропускания</param>
    /// <param name="Gp">Коэффициент передачи в полосе пропускания, по умолчанию 0.891250938 (эквивалент -1 дБ).
    /// Должен быть в диапазоне (0, 1]</param>
    /// <param name="Gs">Коэффициент передачи в полосе заграждения, по умолчанию 0.01 (эквивалент -40 дБ).
    /// Должен быть в диапазоне (0, 1) и меньше Gp</param>
    /// <exception cref="InvalidOperationException">
    /// Выбрасывается если частота пропускания больше либо равна половине частоты дискретизации,
    /// или если Gp &lt;= Gs
    /// </exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Выбрасывается если любой из параметров выходит за допустимые границы
    /// </exception>
    /// <remarks>
    /// Конструктор автоматически:
    /// - Создает спецификацию фильтра на основе входных параметров
    /// - Вычисляет порядок и коэффициенты фильтра
    /// - Инициализирует базовый класс с полученными коэффициентами
    /// 
    /// Типичные значения Gp и Gs:
    /// - Gp = 0.891250938 соответствует -1 дБ затухания в полосе пропускания
    /// - Gs = 0.01 соответствует -40 дБ затухания в полосе заграждения
    /// </remarks>
    /// <example>
    /// <![CDATA[
    /// // Создание фильтра с частотой дискретизации 1000 Гц
    /// double dt = 0.001; // 1000 Гц
    /// var filter = new ButterworthLowPass(
    ///     dt: dt,
    ///     fp: 100,  // граница пропускания 100 Гц
    ///     fs: 200   // граница заграждения 200 Гц
    /// );
    /// 
    /// // Фильтрование сигнала
    /// double[] signal = new double[] { 1.0, 2.0, 3.0, 2.5, 1.0 };
    /// double[] filtered = signal.Select(filter.Process).ToArray();
    /// ]]>
    /// </example>
    public ButterworthLowPass(double dt, double fp, double fs, double Gp = 0.891250938, double Gs = 0.01)
        : this(GetSpecification(dt, fp, fs, Gp, Gs)) { }

    /// <summary>Инициализирует новый низкочастотный фильтр Баттерворта из спецификации</summary>
    /// <param name="Spec">Спецификация фильтра содержащая все требуемые параметры частотной характеристики</param>
    /// <exception cref="ArgumentNullException">Выбрасывается если Spec содержит некорректные значения</exception>
    /// <remarks>
    /// Этот конструктор позволяет повторно использовать спецификацию и создавать несколько экземпляров
    /// фильтра с одинаковыми параметрами. Полезен при параллельной обработке или когда спецификация
    /// вычисляется отдельно от создания фильтра.
    /// </remarks>
    /// <example>
    /// <![CDATA[
    /// // Создание спецификации
    /// var spec = new Specification(dt: 0.001, fp: 100, fs: 200, Gp: 0.891250938, Gs: 0.01);
    /// 
    /// // Создание фильтра из спецификации
    /// var filter = new ButterworthLowPass(spec);
    /// ]]>
    /// </example>
    public ButterworthLowPass(Specification Spec) : this(GetPolynoms(Spec), Spec) { }

    private ButterworthLowPass((double[] A, double[] B) config, Specification Spec) : base(config.B, config.A, Spec) { }
}