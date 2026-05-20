module MyMod
function println(xs...)
    P = get(task_local_storage(), :pix, nothing)
    if P !== nothing
        Base.println(xs..., " (pix ", P, ")")
    else
        Base.println(xs...)
    end
end

function do_print()
    println("Hello world")
    task_local_storage(:pix, 123) do
        println("Processing something")
    end
end
end

MyMod.do_print()
