using System.Globalization;
using System.Reflection.Emit;
using System.Reflection;

namespace MathCore.DSP;

public class Summator
{
    public static Func<int[], int> CreateDelegateInt()
    {
        var method = new DynamicMethod
        (
            name: "Sum",
            returnType: typeof(int),
            parameterTypes: new[]
            {
                typeof(int[]),
            },
            m: Assembly.GetExecutingAssembly().ManifestModule
        );

        method.DefineParameter(1, ParameterAttributes.None, "items");
        //method.DefineParameter(2, ParameterAttributes.Out, "Value");

        var generator = method.GetILGenerator();

        //generator.DeclareLocal(typeof(int).MakePointerType(), pinned: true);

        generator.DeclareLocal(typeof(int));        // v0 - i
        generator.DeclareLocal(typeof(int));        // v1 - count

        generator.Emit(OpCodes.Ldc_I4, 0);
        generator.Emit(OpCodes.Stloc_0);

        generator.Emit(OpCodes.Ldarg_0);
        generator.Emit(OpCodes.Ldlen);
        generator.Emit(OpCodes.Conv_I4);
        generator.Emit(OpCodes.Stloc_1);
        
        generator.Emit(OpCodes.Ldc_I4, 0);          // [0] :: [s]

        var label_check = generator.DefineLabel();
        generator.Emit(OpCodes.Br_S, label_check);  // goto Check

        var label_iteration = generator.DefineLabel();
        generator.MarkLabel(label_iteration); // Begin:

        //                                          // [s]
        generator.Emit(OpCodes.Ldarg_0);            // [items,s]
        generator.Emit(OpCodes.Ldloc_0);            // [i,items,s]
        generator.Emit(OpCodes.Ldelem_I4);          // [items[i],s]
        generator.Emit(OpCodes.Add);                // [s+items[i]]

        generator.Emit(OpCodes.Ldloc_0);            // [i,s]
        generator.Emit(OpCodes.Ldc_I4, 1);          // [1,i,s]
        generator.Emit(OpCodes.Add);                // [i+1,s]
        generator.Emit(OpCodes.Stloc_0);

        generator.MarkLabel(label_check);           // Check:
        generator.Emit(OpCodes.Ldloc_0);            // [i,s]
        generator.Emit(OpCodes.Ldloc_1);            // [count,i,s]
        generator.Emit(OpCodes.Blt_S, label_iteration); // i < count ? [1,s] : [0,s]
        //                                          // [s]

        generator.Emit(OpCodes.Ret);

        var summator = (Func<int[], int>)method.CreateDelegate(typeof(Func<int[], int>));

        return summator;
    }
}