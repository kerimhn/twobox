import { Navbar }          from "@/components/layout/Navbar";
import { Footer }          from "@/components/layout/Footer";
import { HeroSection }     from "@/components/sections/HeroSection";
import { HowItWorks }      from "@/components/sections/HowItWorks";
import { ExplorerSection } from "@/components/sections/ExplorerSection";
import { WhyForny }        from "@/components/sections/WhyForny";
import { ContactSection }  from "@/components/sections/ContactSection";

export default function Home() {
  return (
    <>
      <Navbar />
      <main>
        <HeroSection />
        <HowItWorks />
        <ExplorerSection />
        <WhyForny />
        <ContactSection />
      </main>
      <Footer />
    </>
  );
}
