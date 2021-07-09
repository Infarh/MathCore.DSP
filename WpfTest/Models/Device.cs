namespace WpfTest.Models
{
    internal class Device
    {
        public int Index { get; }
        public string Name { get; }
        public int Channels { get; }

        public Device(int Index, string Name, int Channels)
        {
            this.Index = Index;
            this.Name = Name;
            this.Channels = Channels;
        }
    }
}
