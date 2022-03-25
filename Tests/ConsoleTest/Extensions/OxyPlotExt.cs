﻿using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.ImageSharp;
using OxyPlot.Series;

namespace MathCore.DSP.Extensions;

public static class OxyPlotExt
{
    public static PlotModel Grid(this PlotModel model)
    {
        if (model.Axes.Count > 0)
            foreach (var axe in model.Axes)
            {
                axe.MajorGridlineColor = OxyColors.Gray;
                axe.MinorGridlineColor = OxyColors.LightGray;

                axe.MajorGridlineStyle = LineStyle.Dash;
                axe.MinorGridlineStyle = LineStyle.Dot;
            }
        else
        {
            model.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Bottom,
                MajorGridlineColor = OxyColors.Gray,
                MinorGridlineColor = OxyColors.LightGray,
                MajorGridlineStyle = LineStyle.Dash,
                MinorGridlineStyle = LineStyle.Dot,
            });
            model.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Left,
                MajorGridlineColor = OxyColors.Gray,
                MinorGridlineColor = OxyColors.LightGray,
                MajorGridlineStyle = LineStyle.Dash,
                MinorGridlineStyle = LineStyle.Dot,
            });
        }

        return model;
    }

    public static PlotModel Line(this PlotModel model, IEnumerable<double> Values, double dx, double x0 = 0)
    {
        model.Series.Add(new LineSeries
        {
            ItemsSource = Values.Select((y, i) => new DataPoint(x0 + i * dx, y)),
        });

        return model;
    }

    public static PlotModel Line(this PlotModel model, IEnumerable<double> Values, OxyColor Color, double dx, double x0 = 0)
    {
        model.Series.Add(new LineSeries
        {
            ItemsSource = Values.Select((y, i) => new DataPoint(x0 + i * dx, y)),
            Color = Color,
        });

        return model;
    }

    public static PlotModel SetBackground(this PlotModel model, OxyColor Color)
    {
        model.Background = Color;

        return model;
    }

    public static FileInfo ToPNG(this PlotModel model, string FilePath, int Width = 800, int Height = 600, double Resolution = 96)
    {
        var file = new FileInfo(FilePath);

        var png_exporter = new PngExporter(Width, Height, Resolution);
        using var png = file.Create();
        png_exporter.Export(model, png);

        return file;
    }
}