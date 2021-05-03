using System;
using System.Globalization;
using System.Windows.Data;

namespace WpfTest.Infrastructure.Converters
{
    class Transparent : IValueConverter
    {
        public object Convert(object v, Type tt, object parameter, CultureInfo c) => v;

        public object ConvertBack(object v, Type tt, object parameter, CultureInfo c) => v;
    }
}
