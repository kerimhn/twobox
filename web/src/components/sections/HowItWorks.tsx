import { Search, BarChart2, Hammer } from "lucide-react";

const steps = [
  {
    icon: Search,
    number: "01",
    title: "Vi analyserer boligen din",
    description:
      "Fortell oss om boligen din – alder, størrelse og oppvarmingssystem. Vi vurderer hvilke tiltak som er aktuelle.",
    color: "text-blue-600",
    bg: "bg-blue-50",
  },
  {
    icon: BarChart2,
    number: "02",
    title: "Vi rangerer tiltakene",
    description:
      "Du får en prioritert liste over energitiltak sortert etter lønnsomhet. Vi viser kostnad, besparelse og tilbakebetalingstid.",
    color: "text-green-700",
    bg: "bg-green-50",
  },
  {
    icon: Hammer,
    number: "03",
    title: "Vi kobler deg med byggmester",
    description:
      "Vi setter deg i kontakt med sertifiserte håndverkere som kan gjennomføre tiltakene – og hjelper med Enova-søknad.",
    color: "text-amber-600",
    bg: "bg-amber-50",
  },
];

export function HowItWorks() {
  return (
    <section id="tjenester" className="py-24 bg-white">
      <div className="container mx-auto px-5 md:px-8">
        <div className="text-center mb-16">
          <span className="text-sm font-semibold uppercase tracking-widest text-green-700">Slik fungerer det</span>
          <h2 className="text-4xl font-bold text-slate-900 mt-2 mb-4">
            Fra analyse til ferdig håndverk
          </h2>
          <p className="text-slate-500 max-w-xl mx-auto text-lg">
            Vi gjør det enkelt å ta de riktige valgene for boligen din.
          </p>
        </div>

        <div className="grid md:grid-cols-3 gap-8 max-w-5xl mx-auto">
          {steps.map((step, i) => (
            <div key={i} className="relative">
              {/* Connector line */}
              {i < steps.length - 1 && (
                <div className="hidden md:block absolute top-10 left-[calc(50%+2.5rem)] w-[calc(100%-5rem)] h-px bg-slate-200 z-0" />
              )}

              <div className="relative z-10 flex flex-col items-center text-center p-6">
                <div className={`w-16 h-16 ${step.bg} rounded-2xl flex items-center justify-center mb-5`}>
                  <step.icon size={26} className={step.color} />
                </div>
                <span className="text-xs font-bold text-slate-300 tracking-widest mb-2">{step.number}</span>
                <h3 className="text-lg font-bold text-slate-900 mb-3">{step.title}</h3>
                <p className="text-slate-500 text-sm leading-relaxed">{step.description}</p>
              </div>
            </div>
          ))}
        </div>
      </div>
    </section>
  );
}
