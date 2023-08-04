namespace MathCore.DSP;

public class FilterTest
{
    public static void Run()
    {
        const double f0 = 16e3; // Задаём частоту несущего колебания 16кГц

        // Задаём частоту дискретизации гарантировано превышающую в два раза частоту несущей (т.Котельникова)
        const double fd = 40e3;
        const double dt = 1 / fd; // Период дискретизации

        // Определяем параметры фильтра
        // Сперва аналоговый прототип ФНЧ
        const double fp = 10; // Граничная частота полосы пропускания 10Гц
        const double fs = 40; // Граничная частота полосы подавления 40 Гц
        const double Rp = 1; // Неравноперность АЧХ в полосе пропускания не более 1дБ
        const double Rs = 30; // Уровень подавления 30 дБ в полосе подавления

        // Параметры нужны для определения порядка фильтра:
        // Чем уже полоса между граничной частотой пропускания и заграждения, тем выше порядок
        // Чем меньше величина неоднородности в полосе пропускания и чем глубже подавление в полосе подавления, тем выше порядок

        // Определяем параметры АЧХ фильтра в  единицах
        var Gp = Math.Pow(10.0, -Rp / 20.0);
        var Gs = Math.Pow(10.0, -Rs / 20.0);
        Console.WriteLine("Уровень АЧХ в полосе пропускания не ниже {0:f2}, в полосе заграждения не выше {1:f4}", Gp, Gs);

        // Определяем минимальную и максимальную частоты полосы пропускания фильтра (ширина полосы 1кГц)
        const double f_min = 15.5 * 1e3;
        const double f_max = 16.5 * 1e3;

        // Определяем частоты гармоник слева и справа от несущей (больше по частоте и меньше) для демонстрации эффекта фильтрации
        const double f1 = 1e3;
        const double f2 = 19e3;

        // Создаём фильтр
        //var filter = new ButterworthBandPass(dt, fp, fs, f_min, f_max);
        var filter = new ButterworthBandPass(
            dt,             // Указываем период дискретизации
            fp, fs,         // Указываем параметры ФНЧ-прототипа
            f_min, f_max,   // Указываем желаемую полосу пропускания
            Gp, Gs          // Указываем параметры АЧХ в полосе пропускания и подавления (опционально - если не указать, то будут использованы 1 и 30 дБ соответственно)
        );

        Console.WriteLine("Порядок фильтра {0}", filter.Order);

        // Фильтр умеет фильтровать сам, но можно и вручную.
        // Для этого достаём из фильтра коэффициенты 
        var a = filter.A.ToArray(); // полинома знаменателя
        var b = filter.B.ToArray(); // и полинома числителя

        // Пусть у нас будет сигнал длиной 40 000 отсчётов
        const int samples_count = 40000;

        // Создаём массивы исходных сигналов
        var s0 = new double[samples_count]; // Чистая синусоида на несущей (центральной) частоте 16кГц
        var s1 = new double[samples_count]; // Синусоида частотой ниже несущей (вне полосы пропускания)
        var s2 = new double[samples_count]; // Синусоида частотой выше несущей (вне полосу пропускания)
        var s = new double[samples_count];  // Будет суммарный сигнал

        // Заполняем массивы
        for (var i = 0; i < samples_count; i++)
        {
            var t = i * dt;
            s0[i] = Math.Sin(2 * Math.PI * f0 * t);
            s1[i] = Math.Sin(2 * Math.PI * f1 * t);
            s2[i] = Math.Sin(2 * Math.PI * f2 * t);

            // Так будет рассчитаны значения суммарного сигнала
            s[i] = s0[i] + s1[i] + s2[i];
        }

        var (min_0, max_0, middle_0, delta_0) = GetMinMax(s0);
        var (min_1, max_1, middle_1, delta_1) = GetMinMax(s1);
        var (min_2, max_2, middle_2, delta_2) = GetMinMax(s2);
        var (min_, max_, middle_, delta_) = GetMinMax(s);

        Console.WriteLine("До фильтрации:");
        Console.WriteLine("Амплитуда на несущей частоте      {0}", delta_0);
        Console.WriteLine("Амплитуда на частоте ниже несущей {0}", delta_1);
        Console.WriteLine("Амплитуда на частоте выше несущей {0}", delta_2);
        Console.WriteLine("Амплитуда суммарного сигнала      {0}", delta_);

        // Определяем массивы сигналов на выходе фильтра. Соответственно...
        var y0 = new double[samples_count]; // сигнал на центральной частоте 16 кГц
        var y1 = new double[samples_count]; // сигнал ниже полосы пропускания
        var y2 = new double[samples_count]; // сигнал выше полосы пропускания
        var y = new double[samples_count];  // суммарный сигнал

        // Если мы хотим фильтровать сигналы независимо (как бы по отдельному фильтру на каждый сигнал),
        // то надо создать массив состояния для каждого фильтра
        var filter_state0 = new double[a.Length];
        var filter_state1 = new double[a.Length];
        var filter_state2 = new double[a.Length];
        var filter_state = new double[a.Length];
        // Размер массива должен быть равен порядку фильтра + 1 = длина массива коэффициентов знаменателя (А)

        for (var i = 0; i < samples_count; i++)
        {
            //y0[0] = filter.Process(s0[0]); // Так можно фильтровать отсчёты сигнала самим фильтром. У него внутри есть свой вектор состояния

            // Функция DoubleArrayDSPExtensions.FilterSample реализует логику фильтраиции БИХ фильтра
            // Надо передать в параметрах массив с вектором состояния фильтра, коэффициенты знаменателя и числителя, а также фильтруемое значение сигнала.
            y0[i] = DoubleArrayDSPExtensions.FilterSample(filter_state0, a, b, s0[i]);
            y1[i] = DoubleArrayDSPExtensions.FilterSample(filter_state1, a, b, s1[i]);
            y2[i] = DoubleArrayDSPExtensions.FilterSample(filter_state2, a, b, s2[i]);
            y[i] = DoubleArrayDSPExtensions.FilterSample(filter_state, a, b, s[i]);
        }

        // Отладчиком видно что амплитуда сигнала...
        var (min0, max0, middle0, delta0) = GetMinMax(y0); // на несущей не пострадала
        var (min1, max1, middle1, delta1) = GetMinMax(y1); // вне полосы пропускания была подавлена
        var (min2, max2, middle2, delta2) = GetMinMax(y2);
        var (min, max, middle, delta) = GetMinMax(y);      // суммы равна амплитуде сигнала несущей (не увеличиласть на амплитуду внеполосных сигналов)

        Console.WriteLine();
        Console.WriteLine("После фильтрации:");
        Console.WriteLine("Амплитуда на несущей частоте      {0}", delta0);
        Console.WriteLine("Амплитуда на частоте ниже несущей {0}", delta1);
        Console.WriteLine("Амплитуда на частоте выше несущей {0}", delta2);
        Console.WriteLine("Амплитуда суммарного сигнала      {0}", delta);

        // Спектр также можно посмотреть если интересно.
    }

    private static (double min, double max, double middle, double delta) GetMinMax(double[] array)
    {
        if (array is null) throw new ArgumentNullException(nameof(array));
        if (array.Length == 0) return (double.NaN, double.NaN, double.NaN, double.NaN);
        if (array.Length == 1) return (array[0], array[0], array[0], 0);

        var min = array[0];
        var max = array[0];

        for (var i = 1; i < array.Length; i++)
        {
            var x = array[i];
            if (x < min) min = x;
            if (x > max) max = x;
        }

        return (min, max, (min + max) / 2, max - min);
    }
}